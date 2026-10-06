######################################################################
#   DIANN ➜ limma pipeline – Block × TissueType × Antibody design    #
######################################################################

###############################
## 0. CONFIG ##################
###############################


wd       = r"()" #working directory
pg_file  = "report.pg_matrix.tsv"

CFG <- list(
  wd       = wd,
  pg_file  = pg_file,
  pr_file  = "report.pr_matrix.tsv",
  out_file = "new_egfr_tumor_block_additive_nonorm_noimpute_final.xlsx",
  
  # ---------- filters ----------
  min_peptides     = 2,
  missing_mode     = "global",   # "none" | "global" | "per_arm"
  max_missing_frac = 0.8,
  
  # ---------- processing ----------
  norm_mode   = "none",   # "vsn_global" | "vsn_by_group" | "none"
  impute_mode = "none",           # "MinProb" | "none"
  
  # ---------- modelling ----------
  model_type     = "additive", # "interaction" | "additive" | "two_group"
  block_strategy = "block",       # "block" | "fixed" | "none"
  
  # ---------- hard‑coded column indices ----------
  # Top‑level names = Block / Batch IDs
  # Sub‑list names must contain both the Tissue tag and the Ab tag
  #   e.g. "Healthy_Primary", "Tumor_Isotype"
  # Healthy refers to adjacent non-tumor tissue but is not truly healthy
  # Tumor refers to tumor biopsy tissue
  
  idx = list(
    # ───────── Patient 1 ─────────
    B01 = list(                # Tumor block
      Tumor_Primary  = c( 6, 8,10,12,14,16),
      Tumor_Isotype  = c( 5, 7, 9,11,13,15)
    ),
    B02 = list(                # Adjacent block
      Healthy_Primary = c(18,20,22,24,26,28),
      Healthy_Isotype = c(17,19,21,23,25,27)
    ),
    
    # ───────── Patient 2 ─────────
    B03 = list(
      Tumor_Primary   = c(30,32,34,36,38,40),
      Tumor_Isotype   = c(29,31,33,35,37,39)
    ),
    B04 = list(
      Healthy_Primary = c(42,44,46,48,50,52),
      Healthy_Isotype = c(41,43,45,47,49,51)
    ),
    
    # ───────── Patient 3 ─────────
    B05 = list(
      Healthy_Primary   = c(54,56,58,60,62,64),
      Healthy_Isotype   = c(53,55,57,59,61,63)
    ),
    B06 = list(
      Tumor_Primary = c(66,68,70,72,74,76),
      Tumor_Isotype = c(65,67,69,71,73,75)
    ),
    
    # ───────── Patient 4 ─────────
    B07 = list(
      Healthy_Primary   = c(78,80,82,84,86,88),
      Healthy_Isotype   = c(77,79,81,83,85,87)
    ),
    B08 = list(
      Tumor_Primary = c(90,92,94,96,98,100),
      Tumor_Isotype = c(89,91,93,95,97,99)
    ),
    
    # ───────── Patient 5 ─────────
    B09 = list(
      Healthy_Primary   = c(102,104,106,108,110,112),
      Healthy_Isotype   = c(101,103,105,107,109,111)
    ),
    B10 = list(
      Tumor_Primary = c(114,116,118,120,122,124),
      Tumor_Isotype = c(113,115,117,119,121,123)
    ),
    
    # ───────── Patient 6 ─────────
    B11 = list(
      Tumor_Primary   = c(125,127,129,131,133,135),
      Tumor_Isotype   = c(126,128,130,132,134,136)
    ),
    B12 = list(
      Healthy_Primary = c(137,139,141,143,145,147),
      Healthy_Isotype = c(138,140,142,144,146,148)
    ),
    
    #───────── Patient 7 ─────────
    B13 = list(
      Healthy_Primary   = c(162,164,166,168,170,172),
      Healthy_Isotype   = c(161,163,165,167,169,171)
    ),
    B14 = list(
      Tumor_Primary = c(150,152,154,156,158,160),
      Tumor_Isotype = c(149,151,153,155,157,159)
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
  meta$Ab     <- factor(meta$Ab, levels = c("Primary","Isotype"))
  meta$Block  <- factor(meta$Block)
  
  base_formula <- if(cfg$model_type=="interaction") {
    ~ Tissue * Ab
  } else if(cfg$model_type=="additive") {
    ~ Tissue + Ab
  } else if(cfg$model_type=="two_group") {
    ~ Ab
  } else {
    stop("model_type")
  }
  
  if(cfg$block_strategy == "fixed") {
    rhs <- as.character(base_formula)[2]
    design <- model.matrix(as.formula(paste("~ Block +", rhs)), data = meta)
    eBayes(lmFit(mat, design))
  } else if(cfg$block_strategy == "block") {
    design <- model.matrix(base_formula, data = meta)
    corfit <- duplicateCorrelation(mat, design, block = meta$Block)
    print("block effect")
    print(corfit$consensus)    # e.g. 0.015  → negligible intra‑block correlation
    eBayes(lmFit(mat, design, block = meta$Block,
                 correlation = corfit$consensus))
  } else if(cfg$block_strategy == "none") {
    design <- model.matrix(base_formula, data = meta)
    eBayes(lmFit(mat, design))
  } else {
    stop("block_strategy")
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
  int_col     <- if(!is.na(tissue_col) && !is.na(ab_col)) {
    int_pattern <- paste0("^", tissue_col, ":", ab_col)     # interaction name
    grep(int_pattern, cn, value = TRUE)
  } else {
    character(0)
  }
  
  get_col <- function(mat, colname){
    if(!is.na(colname) && colname %in% colnames(mat)){
      mat[, colname, drop = TRUE]
    } else {
      rep(NA_real_, nrow(mat))
    }
  }
  
  tibble(
    Protein      = rownames(co),
    
    Beta_Tissue  = get_col(co, tissue_col),
    Beta_Ab      = get_col(co, ab_col),
    Beta_Int     = if (length(int_col)) get_col(co, int_col[1]) else NA_real_,
    
    P_Tissue     = get_col(pv, tissue_col),
    P_Ab         = get_col(pv, ab_col),
    P_Int        = if (length(int_col)) get_col(pv, int_col[1]) else NA_real_,
    
    Q_Tissue     = get_col(qv, tissue_col),
    Q_Ab         = get_col(qv, ab_col),
    Q_Int        = if (length(int_col)) get_col(qv, int_col[1]) else NA_real_
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
