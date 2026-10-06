# ==========================================
# Libraries
# ==========================================
library("MSnbase")
library("MSstats")
library("limma")
library("EnhancedVolcano")
library("DEP")
library("dplyr")

# ==========================================
# User Parameters
# ==========================================
working_directory <- "" #file path to DIA-NN output folder


do_imputation <- FALSE
make_files <- TRUE
imputation_filter <- "none"
peptide_count_cutoff <- 0
data_identifier <-'' #common characters in all file names (e.g. initials)
abundance_cutoff <- 0
proteins_of_interest <- c()

proteins_of_interest <- c() #E.g. 'PTPRC'

min_replicates <- 2  # <- Minimum number of non-missing replicates per arm

# Arm definitions
group_1 <- c(6,8,10,12,14,16,18,20,22,24,26,28) #CD3 Neg
group_2 <- c(5,7,9,11,13,15,17,19,21,23,25,27)  #CD3 Pos
group_3 <- c(29,31,33,35,37,39) #CD4 Neg
group_4 <- c(30,32,34,36,38,40) #CD4 Pos
group_5 <- c(41,43,45,47,49,51) #CD8 Neg
group_6 <- c(42,44,46,48,50,52) #CD8 Po
group_7 <- c(53,55,57,59,61,63) #CD20 Neg
group_8 <- c(54,56,58,60,62,64) #CD20 Pos
group_9 <- c(65,67,69,71,73,75) #CD45 Neg
group_10 <- c(66,68,70,72,74,76) #CD45 Pos
group_11 <- c(77,79,81,83,85,87) #MHC Neg
group_12 <- c(78,80,82,84,86,88) #MHC Pos

#H1
group_1H <- c(6,8,10,12,14,16)
group_2H <- c(5,7,9,11,13,15)

pos_column_names <- group_10
neg_column_names <- group_9

pos_group_name <- "" #e.g. 'CD8 Targeted'
neg_group_name <- "" #e.g. 'Isotype'

# Volcano plot cutoffs
PValue_cutoff <- 0.1
logFC_cutoff <- 0.75

data_type <- "DIA"
processing_type <- "separate_arms"  # Options: "none", "uniform", "separate_arms"

# ==========================================
# Directory & File Handling
# ==========================================
if(data_type == "DDA") {
  gene_column <- "Gene"
  protein_column <- "Protein.ID"
  metadata_columns <- 11
  filename <- "%scombined_protein.tsv"
}
if(data_type == "DIA") {
  gene_column <- "Genes"
  protein_column <- "Protein.Group"
  metadata_columns <- 4
  filename <- "%sreport.pg_matrix.tsv"
}

# ==========================================
# Load Data
# ==========================================
# Define data file path
filename <- "report.pg_matrix.tsv"
data_path <- file.path(working_directory, filename)

# Read the data
data <- read.csv(data_path, sep = "\t") 

# Process file to add column names
data_unique <- make_unique(data, gene_column, protein_column)
rownames(data_unique) <- data_unique[, "name"]

names_column <- grep(data_identifier, colnames(data_unique))
data_unique[,names_column][data_unique[,names_column] == 0] <- NA

# For DIA, get peptide counts
if(data_type == "DIA") {
  data_peptides <- read.csv(sprintf("%sreport.pr_matrix.tsv", working_directory), sep = "\t")
  unique_peptides <- distinct(data_peptides)
  unique_peptides <- unique_peptides[order(unique_peptides$Stripped.Sequence),]
  data_peptides_counted <- add_count(unique_peptides, Protein.Group)
  data_peptides_counted <- distinct(data_peptides_counted, Genes, Protein.Group, n)
  data_peptides_counted_unique <- as.data.frame(make_unique(data_peptides_counted, "Genes", "Protein.Group"))
  colnames(data_peptides_counted_unique)[3] <- "Combined.Total.Peptides"
  rownames(data_peptides_counted_unique) <- data_peptides_counted_unique[,"name"]
  data_unique <- merge(data_unique, 
                       data_peptides_counted_unique[3],
                       by = 0)
  data_unique <- data_unique[,-1]
  rownames(data_unique) <- data_unique[,"name"]
}

# Rename arm-containing columns 
colnames(data_unique)[pos_column_names] <- "high_a"
colnames(data_unique)[neg_column_names] <- "low_a"

# ==========================================
# Filter by peptide count
# ==========================================
data_unique_filtered <- subset(data_unique, Combined.Total.Peptides > peptide_count_cutoff)

# ==========================================
# Replicate Filter: min_replicates per arm
# ==========================================
# Count non-missing values
neg_non_missing <- rowSums(!is.na(data_unique_filtered[, neg_column_names]))
pos_non_missing <- rowSums(!is.na(data_unique_filtered[, pos_column_names]))

# Apply filtering
if (identical(min_replicates, "all")) {
  # Must be non-missing in *all* replicates for both groups
  data_unique_filtered <- data_unique_filtered[
    neg_non_missing == length(neg_column_names) &
      pos_non_missing == length(pos_column_names),
  ]
} else if (is.numeric(min_replicates) && min_replicates > 0) {
  # Must meet minimum number of replicates per group
  data_unique_filtered <- data_unique_filtered[
    neg_non_missing >= min_replicates &
      pos_non_missing >= min_replicates,
  ]
} else {
  stop("min_replicates must be either a positive integer or 'all'")
}

# ==========================================
# Processing: Normalization & Imputation
# ==========================================

if (processing_type == "separate_arms") {
  msnset_neg <- readMSnSet2(data_unique_filtered, neg_column_names, "name")
  msnset_pos <- readMSnSet2(data_unique_filtered, pos_column_names, "name")
  
  # VSN normalization per arm
  msnset_neg_norm <- normalise(msnset_neg, "vsn")
  msnset_pos_norm <- normalise(msnset_pos, "vsn")
  
  # Combine normalized MSnSets (no imputation yet)
  combined_preprocessing_postnorm <- MSnbase::combine(msnset_neg_norm, msnset_pos_norm)
  
  # Export normalized, non-imputed report
  data_combined_preprocessing_postnorm <- cbind(
    combined_preprocessing_postnorm@featureData@data[1:metadata_columns],
    combined_preprocessing_postnorm@assayData[["exprs"]]
  )
  data_preprocessed_postnorm_unique <- make_unique(
    data_combined_preprocessing_postnorm, gene_column, protein_column
  )
  
  if (do_imputation) {
    # OPTIONAL: MinProb imputation if you ever want to turn it back on
    msnset_neg_imp <- MSnbase::impute(msnset_neg_norm, "MinProb")
    msnset_pos_imp <- MSnbase::impute(msnset_pos_norm, "MinProb")
    combined_preprocessing <- MSnbase::combine(msnset_neg_imp, msnset_pos_imp)
  } else {
    # *** NO IMPUTATION: use normalized data with NAs ***
    combined_preprocessing <- combined_preprocessing_postnorm
  }
}

if (processing_type == "uniform") {
  msnset <- readMSnSet2(data_unique_filtered, names_column, "name")
  
  # VSN normalization across all samples
  msnset_norm <- normalise(msnset, "vsn")
  
  # Export normalized, non-imputed report
  data_combined_preprocessing_postnorm <- cbind(
    msnset_norm@featureData@data[1:metadata_columns],
    msnset_norm@assayData[["exprs"]]
  )
  data_preprocessed_postnorm_unique <- make_unique(
    data_combined_preprocessing_postnorm, gene_column, protein_column
  )
  
  if (do_imputation) {
    combined_preprocessing <- MSnbase::impute(msnset_norm, "MinProb")
  } else {
    # *** NO IMPUTATION: use vsn-normalized data with NAs ***
    combined_preprocessing <- msnset_norm
  }
}

if (processing_type == "none") {
  # Log2 only, no normalization, no imputation
  data_unique_filtered[, names_column] <- log2(data_unique_filtered[names_column])
  msnset <- readMSnSet2(data_unique_filtered, names_column, "name")
  msnset_norm <- msnset
  
  data_combined_preprocessing_postnorm <- cbind(
    msnset_norm@featureData@data[1:metadata_columns],
    msnset_norm@assayData[["exprs"]]
  )
  data_preprocessed_postnorm_unique <- make_unique(
    data_combined_preprocessing_postnorm, gene_column, protein_column
  )
  
  # Already no imputation here
  combined_preprocessing <- msnset_norm
} 


# ==========================================
# Export Expression Matrix 
# ==========================================
data_combined_preprocessing_exprs <- combined_preprocessing@assayData[["exprs"]]
data_combined_preprocessing_names <- combined_preprocessing@featureData@data[1:metadata_columns]
data_preprocessed_unique <- make_unique(cbind(data_combined_preprocessing_names, data_combined_preprocessing_exprs), gene_column, protein_column)

# ==========================================
# LIMMA Analysis
# ==========================================
data_limma <- data.matrix(
  data_preprocessed_unique[
    c(grep("low_a",  colnames(data_preprocessed_unique)),
      grep("high_a", colnames(data_preprocessed_unique)))
  ]
)
colnames(data_limma) <- 1:ncol(data_limma)
design <- cbind(grp1=1, grp2=c(rep(0,length(neg_column_names)), rep(1,length(pos_column_names))))
limma_output <- eBayes(lmFit(data_limma, design))
limma_output_toptable <- topTable(limma_output, coef=2, n=Inf)
limma_output_toptable$neg_log_p <- -log10(limma_output_toptable$P.Value)
limma_output_toptable <- cbind(row.names(limma_output_toptable), limma_output_toptable)
colnames(limma_output_toptable)[1] <- "name"

full_report <- merge(limma_output_toptable, 
                     data_unique_filtered[grep("imputed|Combined.Total.Peptides", colnames(data_unique_filtered))],
                     by=0)
full_report <- full_report[,-1]
rownames(full_report) <- full_report$name
full_report <- merge(full_report, 
                     data_preprocessed_unique[c(grep("low_a", colnames(data_preprocessed_unique)), grep("high_a", colnames(data_preprocessed_unique)))],
                     by=0)
full_report <- full_report[,-1]
colnames(full_report)[1] <- "name"
rownames(full_report) <- full_report$name
full_report <- merge(full_report, data_unique_filtered[,1:metadata_columns], by=0)
full_report <- full_report[,-1]
rownames(full_report) <- full_report$name
full_report <- full_report[order(full_report$P.Value),]
full_report <- subset(full_report, AveExpr > abundance_cutoff)
filename <- sprintf("DEA_full_report_%s_vs_%s_%s_%s_samples.csv",
                    pos_group_name, neg_group_name, processing_type, min_replicates)
if (make_files) {write.csv(full_report, file.path(working_directory, filename), row.names = FALSE)
}

# ==========================================
# Volcano Plot
# ==========================================
cat("Starting volcano plot...\n")


# Handle proteins of interest safely
if (length(proteins_of_interest) > 0) {
  poi_positions <- which(full_report$name %in% proteins_of_interest)
  full_report$is_poi <- 0
  full_report$is_poi[poi_positions] <- 1
  
  # Reorder so POIs are at the bottom
  full_report_no_poi <- full_report[-poi_positions,]
  full_report_no_poi$is_poi <- 0
  full_report_poi <- full_report[poi_positions,]
  full_report_poi$is_poi <- 1
  full_report <- rbind(full_report_no_poi, full_report_poi)
} else {
  full_report$is_poi <- 0
}

# ==========================================
# Define colors depending on imputation status
# ==========================================

#color names:


if (do_imputation) {
  # --- Original behavior: color by imputation status + POIs ---
  message("Volcano coloring: using imputation-based colors")
  keyvals <- ifelse(
    full_report$name %in% proteins_of_interest, "orange3",
    ifelse(full_report$fully_imputed == 1 & full_report$logFC > logFC_cutoff & full_report$P.Value < PValue_cutoff, "firebrick1",
           ifelse(full_report$partially_imputed == 1 & full_report$logFC > logFC_cutoff & full_report$P.Value < PValue_cutoff, "firebrick3",
                  ifelse(full_report$not_imputed == 1 & full_report$logFC > logFC_cutoff & full_report$P.Value < PValue_cutoff, "darkred", "gray67")))
  )
  keyvals[is.na(keyvals)] <- "black"
  
  unique_colors <- unique(keyvals)
  names(keyvals)[keyvals == "firebrick1"] <- "fully imputed"
  names(keyvals)[keyvals == "firebrick3"] <- "partially imputed"
  names(keyvals)[keyvals == "darkred"]    <- "not imputed"
  if ("orange3" %in% unique_colors) {
    names(keyvals)[keyvals == "orange3"] <- "proteins_of_interest"
  }
  
  # Labels = imputation groups (+ POIs if present)
  if (length(proteins_of_interest) > 0) {
    selectLab <- full_report$name[
      which(names(keyvals) %in% c(
        "fully imputed", "partially imputed", "not imputed", "proteins_of_interest"
      ))
    ]
  } else {
    selectLab <- full_report$name[
      which(names(keyvals) %in% c(
        "fully imputed", "partially imputed", "not imputed"
      ))
    ]
  }
  
} else {
  message("Volcano coloring: using NON-IMPUTED sign-based colors")
  
  # Make sure cutoff is positive & symmetric
  if (logFC_cutoff < 0) {
    logFC_cutoff <- abs(logFC_cutoff)
    message("logFC_cutoff was negative; using abs(logFC_cutoff) = ", logFC_cutoff)
  }
  
  # Significant positive / negative based on your thresholds
  sig_pos <- full_report$logFC >  logFC_cutoff  & full_report$P.Value < PValue_cutoff
  sig_neg <- full_report$logFC < -logFC_cutoff  & full_report$P.Value < PValue_cutoff
  
  message("Significant positive proteins: ", sum(sig_pos, na.rm = TRUE))
  message("Significant negative proteins: ", sum(sig_neg, na.rm = TRUE))
  
  
  # Start everything as non-significant
  keyvals <- rep("gray67", nrow(full_report))
  
  #assign colors to each group
  if (pos_group_name == "OvCan") {
    keyvals[sig_pos] <- "deepskyblue"
  } else if (pos_group_name == "Leiomyoma") {
    keyvals[sig_pos] <- "palevioletred1"
  } else if (pos_group_name == "EMCan") {
    keyvals[sig_pos] <- "mediumseagreen"
  } else if (pos_group_name == "AdCyst") {
    keyvals[sig_pos] <- "mediumpurple1"
  }
  
  if (neg_group_name == "OvCan") {
    keyvals[sig_neg] <- "deepskyblue"
  } else if (neg_group_name == "Leiomyoma") {
    keyvals[sig_neg] <- "palevioletred1"
  } else if (neg_group_name == "EMCan") {
    keyvals[sig_neg] <- "mediumseagreen"
  } else if (neg_group_name == "AdCyst") {
    keyvals[sig_neg] <- "mediumpurple1"
  }
  
  # Overwrite with POI color last, so POIs always orange
  if (length(proteins_of_interest) > 0) {
    poi_positions <- which(full_report$name %in% proteins_of_interest)
    keyvals[poi_positions] <- "black"
  }
  
  # IMPORTANT: names = colors themselves, so EnhancedVolcano can map them cleanly
  names(keyvals) <- keyvals
  
  # Label only POIs (you can change this)
  if (length(proteins_of_interest) > 0) {
    selectLab <- full_report$name[full_report$name %in% proteins_of_interest]
  } else {
    selectLab <- character(0)
  }
}

# ==========================================
# Axis limits
# ==========================================
buffer_x <- 0.05 * diff(range(full_report$logFC, na.rm = TRUE))
buffer_y <- 0.05 * diff(range(full_report$neg_log_p, na.rm = TRUE))
plot_min_logfc <- -4
plot_max_logfc <- 4
plot_min_logp  <- 0
plot_max_logp  <- max(full_report$neg_log_p, na.rm = TRUE) + buffer_y

# ==========================================
# Volcano Plot (Label significant proteins)
# ==========================================
cat("Starting volcano plot...\n")

# Determine significance: pass both logFC AND p-value cutoffs
sig_pos <- full_report$logFC >  logFC_cutoff  & full_report$P.Value < PValue_cutoff
sig_neg <- full_report$logFC < -logFC_cutoff  & full_report$P.Value < PValue_cutoff
sig_any <- sig_pos | sig_neg

# Handle proteins of interest safely
full_report$is_poi <- 0
if (length(proteins_of_interest) > 0) {
  poi_positions <- which(full_report$name %in% proteins_of_interest)
  full_report$is_poi[poi_positions] <- 1
}

# Start everything as non-significant
keyvals <- rep("gray67", nrow(full_report))

# Assign colors to significant proteins
keyvals[sig_pos] <- "deepskyblue"  # positive significant
keyvals[sig_neg] <- "red"          # negative significant

# Overwrite POI color last
if (length(proteins_of_interest) > 0) {
  keyvals[full_report$is_poi == 1] <- "black"
}

# IMPORTANT: names = colors themselves for EnhancedVolcano
names(keyvals) <- keyvals

# Label proteins: all significant + POIs
selectLab <- full_report$name[sig_any | full_report$is_poi == 1]

# ==========================================
# Axis limits
# ==========================================
buffer_x <- 0.05 * diff(range(full_report$logFC, na.rm = TRUE))
buffer_y <- 0.05 * diff(range(full_report$neg_log_p, na.rm = TRUE))
plot_min_logfc <- -5
plot_max_logfc <- 5
plot_min_logp  <- 0
plot_max_logp  <- max(full_report$neg_log_p, na.rm = TRUE) + buffer_y

# ==========================================
# Plot volcano
# ==========================================
volcano_plot <- EnhancedVolcano(
  full_report,
  lab = full_report$name,
  x = 'logFC',
  y = 'P.Value',
  xlab = paste("log2 (",pos_group_name,"/",neg_group_name,")"),
  ylab = "-log10(p-value)",
  xlim = c(plot_min_logfc, plot_max_logfc),
  ylim = c(plot_min_logp, plot_max_logp),
  labSize = 4,
  pointSize = ifelse(full_report$is_poi > 0, 3, 3),
  colAlpha = 1,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  arrowheads = FALSE,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  subtitle = NULL,
  legendPosition = "none",
  title = NULL,
  cutoffLineType = "blank",
  colCustom = keyvals,
  selectLab = selectLab
)

volcano_plot <- volcano_plot +
  geom_hline(
    yintercept = -log10(PValue_cutoff),
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  )

print(volcano_plot)

#plot name
plot_name <- sprintf("RPlot_%s_vs_%s_%s_%s_samples.tif",
                     pos_group_name, neg_group_name, processing_type,min_replicates)

# Export volcano as TIF
if (make_files) {tiff(
  file.path(working_directory, plot_name),
  width = 5500,
  height = 4000,
  res = 600,
  compression = "lzw"
)


print(volcano_plot)
dev.off()
}
