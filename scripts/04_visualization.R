# ============================================================
# 04 — VISUALIZATION
# Dynamic region (ITS1, ITS2, ITS) and taxonomic rank selection
# Works with dnabarcoder '.classified' Full classification column
# ============================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 1. Parse command-line arguments
# ------------------------------------------------------------

cli_args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(args, flag, default_val) {
  idx <- grep(paste0("^", flag, "="), args)
  if (length(idx) > 0) {
    return(sub(paste0("^", flag, "="), "", args[idx[1]]))
  }
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) {
    return(args[idx + 1])
  }
  return(default_val)
}

REGION <- tolower(parse_arg(cli_args, "--region", parse_arg(cli_args, "-r", "its2")))
RANK   <- tolower(parse_arg(cli_args, "--rank", parse_arg(cli_args, "-k", "genus")))
REMOVE_UNCLASSIFIED <- "--remove-unidentified" %in% cli_args

valid_regions <- c("its1", "its2", "its")
if (!REGION %in% valid_regions) {
  stop("Invalid --region. Must be one of: ", paste(valid_regions, collapse = ", "))
}

valid_ranks <- c("phylum", "class", "order", "family", "genus", "species")
if (!RANK %in% valid_ranks) {
  stop("Invalid --rank. Must be one of: ", paste(valid_ranks, collapse = ", "))
}

# Prefix mapping for standard UNITE/dnabarcoder taxonomy strings
rank_prefixes <- c(
  phylum  = "p__",
  class   = "c__",
  order   = "o__",
  family  = "f__",
  genus   = "g__",
  species = "s__"
)
target_prefix <- rank_prefixes[RANK]

# Find project folder
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
  PROJECT_ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
} else {
  PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
PLOTS_DIR <- file.path(RESULTS_DIR, "plots")
DNABARCODER_DIR <- file.path(RESULTS_DIR, "dnabarcoder")

dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("VISUALIZATION SETUP\n")
cat("========================================\n")
cat("Project root: ", PROJECT_ROOT, "\n")
cat("Region:       ", toupper(REGION), "\n")
cat("Taxon Rank:   ", RANK, "\n\n")

# ------------------------------------------------------------
# 2. Load packages
# ------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is not installed.")
if (!requireNamespace("vegan", quietly = TRUE)) stop("vegan is not installed.")

library(ggplot2)
library(vegan)

# ------------------------------------------------------------
# 3. Read DADA2 ASV table
# ------------------------------------------------------------

seqtab_file <- file.path(RESULTS_DIR, "seqtab_nochim.rds")

if (!file.exists(seqtab_file)) {
  stop("Could not find:\n", seqtab_file, "\nRun the DADA2 step first.")
}

seqtab <- readRDS(seqtab_file)
if (!is.matrix(seqtab)) seqtab <- as.matrix(seqtab)

cat("ASV table:\n")
cat("  Samples: ", nrow(seqtab), "\n")
cat("  ASVs:    ", ncol(seqtab), "\n\n")

# ------------------------------------------------------------
# 4. Read dnabarcoder taxonomy matching selected region
# ------------------------------------------------------------

taxonomy_filename <- paste0("asvs.unite2025", REGION, "_BLAST.classified")
taxonomy_file <- file.path(DNABARCODER_DIR, taxonomy_filename)

if (!file.exists(taxonomy_file)) {
  stop("Could not find taxonomy file:\n", taxonomy_file, 
       "\nMake sure you ran 03_taxonomy.py with '--region ", REGION, "'.")
}

cat("Taxonomy file:\n  ", taxonomy_file, "\n\n")

taxonomy <- read.delim(
  taxonomy_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

# ------------------------------------------------------------
# 5. Extract requested rank from dnabarcoder output
# ------------------------------------------------------------

if (!"ID" %in% names(taxonomy)) {
  stop("Could not find 'ID' column in taxonomy file.")
}

# Determine taxonomy text column (Full classification, Prediction, or direct rank column)
tax_col <- NULL
if (RANK %in% tolower(names(taxonomy))) {
  tax_col <- names(taxonomy)[tolower(names(taxonomy)) == RANK][1]
} else if ("Full classification" %in% names(taxonomy)) {
  tax_col <- "Full classification"
} else if ("Prediction" %in% names(taxonomy)) {
  tax_col <- "Prediction"
} else {
  stop("Could not find classification details in taxonomy file.")
}

taxonomy_id <- as.character(taxonomy[["ID"]])
raw_tax_values <- as.character(taxonomy[[tax_col]])
asv_names <- colnames(seqtab)

# Function to parse specific rank from standard strings (e.g. k__Fungi;p__Ascomycota;g__Aspergillus)
extract_rank <- function(tax_string, prefix) {
  if (is.na(tax_string) || tax_string == "") return("Unclassified")
  
  # If string directly contains the rank prefix (e.g., g__Aspergillus)
  pattern <- paste0(".*?", prefix, "([^;]+).*")
  if (grepl(pattern, tax_string, ignore.case = TRUE)) {
    val <- sub(pattern, "\\1", tax_string, ignore.case = TRUE)
    return(trimws(val))
  }
  
  # If column was already named as the rank directly without prefix
  if (!grepl(";", tax_string) && !grepl("__", tax_string)) {
    return(trimws(tax_string))
  }
  
  return("Unclassified")
}

taxon_extracted <- sapply(raw_tax_values, extract_rank, prefix = target_prefix, USE.NAMES = FALSE)

# ------------------------------------------------------------
# 6. Match ASVs by number
# ------------------------------------------------------------

asv_number <- suppressWarnings(
  as.integer(sub("^ASV", "", taxonomy_id, ignore.case = TRUE))
)

if (all(!is.na(asv_number)) && all(asv_number >= 1) && all(asv_number <= length(asv_names))) {
  match_index <- asv_number
} else {
  stop("Could not match taxonomy IDs to DADA2 ASVs. Expected IDs such as ASV1, ASV2...")
}

# ------------------------------------------------------------
# 7. Assign taxon to every ASV and clean labels
# ------------------------------------------------------------

taxon_for_asv <- rep("Unclassified", length(asv_names))

valid <- (!is.na(match_index) & match_index >= 1 & match_index <= length(taxon_extracted))
taxon_for_asv[valid] <- taxon_extracted[match_index[valid]]

taxon_for_asv <- trimws(taxon_for_asv)

# Strip any residual prefixes (e.g., g__, f__)
taxon_for_asv <- sub("^[pcofgs]__", "", taxon_for_asv, ignore.case = TRUE)

taxon_for_asv[
  is.na(taxon_for_asv) |
  taxon_for_asv == "" |
  grepl("^unidentified", taxon_for_asv, ignore.case = TRUE) |
  tolower(taxon_for_asv) %in% c(
    "na", "nan", "none", "unknown", "unclassified", 
    "uncultured", "no_match", "incertae_sedis", "incertae sedis"
  )
] <- "Unclassified"

cat("ASVs assigned to ", RANK, ": ", 
    sum(taxon_for_asv != "Unclassified"), " / ", length(taxon_for_asv), "\n\n", sep = "")

# ------------------------------------------------------------
# 8. Combine ASVs into taxonomic rank
# ------------------------------------------------------------

rank_levels <- unique(taxon_for_asv)
rank_table <- matrix(0, nrow = nrow(seqtab), ncol = length(rank_levels))
rownames(rank_table) <- rownames(seqtab)
colnames(rank_table) <- rank_levels

for (i in seq_along(rank_levels)) {
  selected <- which(taxon_for_asv == rank_levels[i])
  if (length(selected) == 1) {
    rank_table[, i] <- seqtab[, selected]
  } else {
    rank_table[, i] <- rowSums(seqtab[, selected, drop = FALSE])
  }
}

rank_table <- as.data.frame(rank_table, check.names = FALSE)

if (REMOVE_UNCLASSIFIED) {
  drop_cols <- tolower(names(rank_table)) %in% c("unidentified", "unclassified")
  if (any(drop_cols)) {
    cat("Removing unidentified/unclassified:", sum(rank_table[, drop_cols, drop = FALSE]), "reads\n\n")
    rank_table <- rank_table[, !drop_cols, drop = FALSE]
  }
  if (ncol(rank_table) == 0) stop("No classified taxa left.")
}

write.table(
  rank_table,
  file = file.path(RESULTS_DIR, paste0(RANK, "_abundance_table_", REGION, ".tsv")),
  sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
)

# ------------------------------------------------------------
# 9. Total Abundance Plot
# ------------------------------------------------------------

total_abundance <- colSums(rank_table, na.rm = TRUE)
taxon_total <- data.frame(
  Taxon = names(total_abundance),
  Abundance = as.numeric(total_abundance),
  stringsAsFactors = FALSE
)
taxon_total <- taxon_total[order(-taxon_total$Abundance), ]

TOP_N <- 25
plot_total <- head(taxon_total, min(TOP_N, nrow(taxon_total)))
plot_total$Taxon <- factor(plot_total$Taxon, levels = plot_total$Taxon)

p_total <- ggplot(plot_total, aes(x = Taxon, y = Abundance)) +
  geom_col(width = 0.75) +
  labs(
    title = paste0("Fungal ", RANK, " abundance (", toupper(REGION), ")"),
    subtitle = paste0("Top ", nrow(plot_total), " ", RANK, " taxa across all samples"),
    x = RANK,
    y = "Abundance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(PLOTS_DIR, paste0(RANK, "_total_abundance_", REGION, ".png")),
  plot = p_total, width = 12, height = 7, dpi = 300
)

# ------------------------------------------------------------
# 10. Relative Abundance Plot
# ------------------------------------------------------------

sample_totals <- rowSums(rank_table)
relative_table <- rank_table

for (i in seq_len(nrow(relative_table))) {
  if (sample_totals[i] > 0) {
    relative_table[i, ] <- relative_table[i, ] / sample_totals[i]
  }
}

relative_long <- do.call(rbind, lapply(seq_len(nrow(relative_table)), function(i) {
  data.frame(
    Sample = rownames(relative_table)[i],
    Taxon = names(relative_table),
    RelativeAbundance = as.numeric(relative_table[i, ]),
    stringsAsFactors = FALSE
  )
}))

top_taxa <- head(taxon_total$Taxon, min(TOP_N, nrow(taxon_total)))
relative_long$Taxon <- ifelse(relative_long$Taxon %in% top_taxa, relative_long$Taxon, "Other")

relative_long <- aggregate(RelativeAbundance ~ Sample + Taxon, data = relative_long, FUN = sum)

p_relative <- ggplot(relative_long, aes(x = Sample, y = RelativeAbundance, fill = Taxon)) +
  geom_col(width = 0.9) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  labs(
    title = paste0("Fungal ", RANK, " relative abundance (", toupper(REGION), ")"),
    subtitle = "Top taxa shown separately",
    x = "Sample",
    y = "Relative abundance",
    fill = RANK
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(PLOTS_DIR, paste0(RANK, "_relative_abundance_", REGION, ".png")),
  plot = p_relative, width = 14, height = 8, dpi = 300
)

# ------------------------------------------------------------
# 11. Bray-Curtis PCoA + Clustering
# ------------------------------------------------------------

sample_totals_pcoa <- rowSums(rank_table, na.rm = TRUE)
keep_samples <- sample_totals_pcoa > 0
rank_table_pcoa <- rank_table[keep_samples, , drop = FALSE]

if (nrow(rank_table_pcoa) >= 3) {
  bray <- vegdist(rank_table_pcoa, method = "bray")
  pcoa <- cmdscale(bray, k = 2, eig = TRUE, add = TRUE)

  pcoa_points <- data.frame(
    Sample = rownames(rank_table_pcoa),
    PCoA1 = pcoa$points[, 1],
    PCoA2 = pcoa$points[, 2],
    stringsAsFactors = FALSE
  )

  N_CLUSTERS <- min(3, nrow(rank_table_pcoa) - 1)
  set.seed(123)
  cluster_result <- kmeans(pcoa_points[, c("PCoA1", "PCoA2")], centers = N_CLUSTERS)
  pcoa_points$Cluster <- factor(cluster_result$cluster)

  positive_eig <- pcoa$eig[pcoa$eig > 0]
  percent1 <- if (length(positive_eig) >= 1) round(100 * positive_eig[1] / sum(positive_eig), 1) else NA
  percent2 <- if (length(positive_eig) >= 2) round(100 * positive_eig[2] / sum(positive_eig), 1) else NA

  p_pcoa <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = Cluster)) +
    geom_point(size = 3) +
    labs(
      title = paste0("Bray-Curtis PCoA - ", RANK, " level (", toupper(REGION), ")"),
      subtitle = paste0("Samples clustered into ", N_CLUSTERS, " groups"),
      x = paste0("PCoA1", ifelse(is.na(percent1), "", paste0(" (", percent1, "%)"))),
      y = paste0("PCoA2", ifelse(is.na(percent2), "", paste0(" (", percent2, "%)"))),
      color = "Cluster"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

  ggsave(
    filename = file.path(PLOTS_DIR, paste0("bray_curtis_pcoa_", RANK, "_", REGION, ".png")),
    plot = p_pcoa, width = 8, height = 6, dpi = 300
  )
}

# ------------------------------------------------------------
# 12. Export Taxonomy Mapping
# ------------------------------------------------------------

taxonomy_export <- data.frame(
  ASV = asv_names,
  Taxon = taxon_for_asv,
  stringsAsFactors = FALSE
)
names(taxonomy_export)[2] <- RANK

write.table(
  taxonomy_export,
  file = file.path(RESULTS_DIR, paste0("taxonomy_table_", RANK, "_", REGION, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("========================================\n")
cat("VISUALIZATION COMPLETE\n")
cat("========================================\n")