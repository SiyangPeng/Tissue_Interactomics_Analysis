######################################################################
#   DIANN ➜ limma pipeline – Block × TissueType × Antibody design    #
######################################################################

###############################
## 0. CONFIG ##################
###############################


wd       = "" #Working directory
pg_file  = "report.pg_matrix.tsv"

CFG <- list(
  wd       = wd,
  pg_file  = pg_file,
  pr_file  = "report.pr_matrix.tsv",
  out_file = "interact_test_nonorm.xlsx",
  
  # ---------- filters ----------
  min_peptides     = 1,
  missing_mode     = "per_arm",   # "none" | "global" | "per_arm"
  max_missing_frac = 0.8,
  
  # ---------- processing ----------
  norm_mode   = "none",   # "vsn_global" | "vsn_by_group" | "none"
  impute_mode = "none",           # "MinProb" | "none"
  
  # ---------- modelling ----------
  model_type     = "additive", # "interaction" | "additive"
  
  # ---------- hard‑coded column indices ----------
  # Top‑level names = Block / Batch IDs
  # Sub‑list names must contain both the Tissue tag and the Ab tag
  #   e.g. "Healthy_Primary", "Tumor_Isotype"
  
  idx = list(
    B01 = list(
      Synapse_Isotype = c(5,7,9,11,13,15,17,19,21,23,25,27),
      Synapse_Primary = c(6,8,10,12,14,16,18,20,22,24,26,28),
      Nonsynapse_Isotype  = c(29,31,33,35,37,39,41,43,45,47,49,51),
      Nonsynapse_Primary  = c(30,32,34,36,38,40,42,44,46,48,50,52)
    )
  )
)

###############################
## 1. LIBRARIES ###############
###############################

suppressPackageStartupMessages({
  library(MSnbase); library(limma); library(dplyr); library(tidyr); library(DEP)
  library(writexl)
})

say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n")

###############################
## 2. SAMPLE METADATA #########
###############################

parse_subname <- function(x){
  # split “Healthy_Primary” → Tissue = Healthy, Ab = Primary
  parts <- strsplit(x, "_", fixed = TRUE)[[1]]
  if(length(parts) != 2) stop("sub‑list names must be Tissue_Ab")
  list(Tissue = parts[1], Ab = parts[2])
}

build_sample_info <- function(idx){
  
  rows <- lapply(names(idx), function(block_id){
    
    lapply(names(idx[[block_id]]), function(sub){
      
      par <- parse_subname(sub)
      cols <- idx[[block_id]][[sub]]
      
      data.frame(
        Col   = cols,
        Block = block_id,
        Tissue = par$Tissue,
        Ab     = par$Ab,
        stringsAsFactors = FALSE
      )
    }) |>
      bind_rows()
  }) |>
    bind_rows() |>
    arrange(Col)
  
  rows$Sample <- with(rows, paste(Block, Tissue, Ab,
                                  sprintf("R%02d", ave(Col, Block, Tissue, Ab, FUN = seq_along)), sep = "_"))
  
  rows
}

###############################
## 3. LOAD & MERGE REPORTS ####
###############################

load_proteins <- function(cfg){
  pg <- read.delim(file.path(cfg$wd, cfg$pg_file), check.names = FALSE)
  pr <- read.delim(file.path(cfg$wd, cfg$pr_file), check.names = FALSE) |>
    distinct(Stripped.Sequence, Protein.Group)
  
  pep_counts <- pr |>
    add_count(Protein.Group, name = "Peptides") |>
    distinct(Protein.Group, Peptides)
  
  pg_u   <- make_unique(pg, "Genes", "Protein.Group")
  merged <- left_join(pg_u, pep_counts, by = "Protein.Group")
  rownames(merged) <- merged$name
  merged
}

###############################
## 4. FILTER / NORMALISE ######
###############################

missing_filter <- function(mat, meta, mode, thr){
  keep <- switch(
    mode,
    
    none   = rep(TRUE, nrow(mat)),
    
    global = rowMeans(is.na(mat)) <= thr,
    
    per_arm = {
      # one key per physical block * tissue * antibody
      keys  <- paste(meta$Block, meta$Tissue, meta$Ab, sep = "_")
      arms  <- split(seq_len(ncol(mat)), keys)
      
      sapply(seq_len(nrow(mat)), function(i)
        all(vapply(arms,
                   function(cols) mean(is.na(mat[i, cols])) <= thr,
                   FUN.VALUE = TRUE)))
    },
    
    stop("missing_mode")
  )
  
  say("Missing‑value filter:", sum(keep), "/", nrow(mat), "proteins kept")
  mat[keep, ]
}

vsn_by_group <- function(msn, meta){
  # one key per physical block AND antibody arm
  grp_key <- paste(meta$Block, meta$Tissue, meta$Ab, sep = "_")
  ex <- exprs(msn)
  
  for (key in unique(grp_key)) {
    cols <- which(grp_key == key)
    ex[, cols] <- exprs(normalise(msn[, cols], "vsn"))
  }
  
  exprs(msn) <- ex
  msn
}

normalise_ms <- function(msn, meta, mode){
  switch(mode,
         none         = msn,
         vsn_global   = normalise(msn, "vsn"),
         vsn_by_group = vsn_by_group(msn, meta),
         stop("norm_mode"))
}

impute_ms <- function(msn, mode){
  switch(mode,
         none    = msn,
         MinProb = MSnbase::impute(msn, "MinProb"),
         stop("impute_mode"))
}

###############################
## 5. LIMMA MODEL #############
###############################

fit_model <- function(mat, meta, cfg){
  
  meta$Tissue <- factor(meta$Tissue)
  meta$Ab     <- factor(meta$Ab, levels = c("Isotype","Primary"))
  meta$Block  <- factor(meta$Block)
  
  design <- if (cfg$model_type == "interaction")
    model.matrix(~ Tissue * Ab, data = meta)
  else
    model.matrix(~ Tissue + Ab, data = meta)
  
  if (nlevels(meta$Block) < 2) {
    ## ---- single block: no random effect --------------------------------
    fit <- lmFit(mat, design)
    eBayes(fit)
  } else {
    ## ---- usual duplicateCorrelation path -------------------------------
    corfit <- duplicateCorrelation(mat, design, block = meta$Block)
    print("block effect")
    print(corfit$consensus)    # e.g. 0.015  → negligible intra‑block correlation
    
    eBayes(lmFit(mat, design,
                 block = meta$Block,
                 correlation = corfit$consensus))
  }
}

tidy_stats <- function(fit){
  
  co <- fit$coefficients
  pv <- fit$p.value
  qv <- apply(pv, 2, p.adjust, "BH")
  cn <- colnames(co)
  
  ## ---- grab column names dynamically ---------------------------------
  tissue_col  <- grep("^Tissue",  cn, value = TRUE)[1]      # first tissue level
  ab_col      <- grep("^Ab",      cn, value = TRUE)[1]      # AbIsotype
  int_pattern <- paste0("^", tissue_col, ":", ab_col)       # interaction name
  int_col     <- grep(int_pattern, cn, value = TRUE)
  
  tibble(
    Protein      = rownames(co),
    
    Beta_Tissue  = co[, tissue_col, drop = TRUE],
    Beta_Ab      = co[, ab_col,     drop = TRUE],
    Beta_Int     = if (length(int_col))
      co[, int_col,  drop = TRUE] else NA_real_,
    
    P_Tissue     = pv[, tissue_col],
    P_Ab         = pv[, ab_col],
    P_Int        = if (length(int_col)) pv[, int_col] else NA_real_,
    
    Q_Tissue     = qv[, tissue_col],
    Q_Ab         = qv[, ab_col],
    Q_Int        = if (length(int_col)) qv[, int_col] else NA_real_
  )
}


`%||%` <- function(x, y) if(is.null(x)) y else x   # helper

###############################
## 6. MAIN ####################
###############################

main <- function(cfg = CFG){
  
  setwd(cfg$wd)
  
  meta <- build_sample_info(cfg$idx)
  rownames(meta) <- meta$Sample
  
  prot <- load_proteins(cfg) |>
    filter(is.na(Peptides) | Peptides >= cfg$min_peptides)
  
  mtx <- as.matrix(prot[, meta$Col]); colnames(mtx) <- meta$Sample
  mtx <- missing_filter(mtx, meta, cfg$missing_mode, cfg$max_missing_frac)
  mtx_pre <- mtx
  
  msn <- MSnSet(exprs = mtx, pData = meta)
  msn <- normalise_ms(msn, meta, cfg$norm_mode)
  if(cfg$norm_mode == "none") exprs(msn) <- log2(exprs(msn))
  msn <- impute_ms(msn, cfg$impute_mode)
  
  fit <- fit_model(exprs(msn), meta, cfg)
  res <- tidy_stats(fit) |>
    mutate(across(starts_with("P_"), ~ -log10(.x),
                  .names = "mLog10_{col}"))
  
  if(cfg$impute_mode == "none"){
    res$ImputeStatus <- "not_imputed"
  } else {
    miss_prop <- rowMeans(is.na(mtx_pre))
    res$ImputeStatus <- ifelse(miss_prop == 1, "fully_imputed",
                               ifelse(miss_prop == 0, "not_needed", "partially_imputed"))
  }
  
  ann_cols <- intersect(c("Description","Genes","Peptides"), names(prot))
  res <- res |>
    left_join(prot[, c("name", ann_cols)], by = c("Protein" = "name")) |>
    left_join(tibble::rownames_to_column(as.data.frame(exprs(msn)), "Protein"),
              by = "Protein")
  
  say("write →", cfg$out_file)
  write_xlsx(list(stats = res), cfg$out_file)
  say("done.")
}

data <- read.csv(file.path(wd, pg_file), sep="\t")
main()
