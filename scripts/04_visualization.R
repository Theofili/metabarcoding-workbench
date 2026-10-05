# ============================================================
# 04 — VISUALIZATION
# Fungal its2 metabarcoding
# ============================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 1. Find project folder
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])

  PROJECT_ROOT <- normalizePath(
    file.path(dirname(script_path), ".."),
    winslash = "/",
    mustWork = TRUE
  )
} else {
  PROJECT_ROOT <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )
}

RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
PLOTS_DIR <- file.path(RESULTS_DIR, "plots")
DNABARCODER_DIR <- file.path(
  RESULTS_DIR,
  "dnabarcoder"
)

dir.create(
  PLOTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

cli_args <- commandArgs(trailingOnly = TRUE)
REMOVE_UNCLASSIFIED <- "--remove-unidentified" %in% cli_args

cat("Project root:\n")
cat(PROJECT_ROOT, "\n\n")


# ------------------------------------------------------------
# 2. Load packages
# ------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is not installed.")
}

if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("vegan is not installed.")
}

library(ggplot2)
library(vegan)


# ------------------------------------------------------------
# 3. Read DADA2 ASV table
# ------------------------------------------------------------

seqtab_file <- file.path(
  RESULTS_DIR,
  "seqtab_nochim.rds"
)

if (!file.exists(seqtab_file)) {
  stop(
    "Could not find:\n",
    seqtab_file,
    "\nRun the DADA2 step first."
  )
}

seqtab <- readRDS(seqtab_file)

if (!is.matrix(seqtab)) {
  seqtab <- as.matrix(seqtab)
}

cat("ASV table:\n")
cat("  Samples: ", nrow(seqtab), "\n", sep = "")
cat("  ASVs: ", ncol(seqtab), "\n\n", sep = "")


# ------------------------------------------------------------
# 4. Read dnabarcoder taxonomy
# ------------------------------------------------------------

taxonomy_file <- file.path(
  DNABARCODER_DIR,
  "asvs.unite2025its2_BLAST.classification"
)

if (!file.exists(taxonomy_file)) {
  stop(
    "Could not find taxonomy file:\n",
    taxonomy_file
  )
}

cat("Taxonomy file:\n")
cat("  ", taxonomy_file, "\n\n", sep = "")

taxonomy <- read.delim(
  taxonomy_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

cat("Taxonomy columns:\n")
cat("  ", paste(names(taxonomy), collapse = ", "), "\n\n",
    sep = "")


# ------------------------------------------------------------
# 5. Find ID and class columns
# ------------------------------------------------------------

if ("ID" %in% names(taxonomy)) {
  id_col <- "ID"
} else {
  stop("Could not find ID column.")
}

if ("genus" %in% names(taxonomy)) {
  genus_col <- "genus"
} else if ("genus" %in% names(taxonomy)) {
  genus_col <- "genus"
} else {
  stop("Could not find genus column.")
}

taxonomy_id <- as.character(
  taxonomy[[id_col]]
)

genus <- as.character(
  taxonomy[[genus_col]]
)

asv_names <- colnames(seqtab)

cat("Taxonomy mapping:\n")
cat("  ID column: ", id_col, "\n", sep = "")
cat("  genus column: ", genus_col, "\n", sep = "")


# ------------------------------------------------------------
# 6. Match ASV1, ASV2, ASV3... to DADA2 ASVs
# ------------------------------------------------------------
#
# dnabarcoder produced IDs such as:
#
#   ASV1
#   ASV2
#   ASV3
#
# DADA2 stores the DNA sequence itself as the ASV name.
#
# Therefore:
#
#   ASV1  -> first DADA2 ASV
#   ASV2  -> second DADA2 ASV
#   ASV3  -> third DADA2 ASV
#
# etc.
# ------------------------------------------------------------

asv_number <- suppressWarnings(
  as.integer(
    sub(
      "^ASV",
      "",
      taxonomy_id,
      ignore.case = TRUE
    )
  )
)

if (
  all(!is.na(asv_number)) &&
  all(asv_number >= 1) &&
  all(asv_number <= length(asv_names))
) {

  match_index <- asv_number

  cat(
    "  IDs are ASV1, ASV2, ASV3, ...\n"
  )

  cat(
    "  Mapping taxonomy IDs to DADA2 ASVs by ASV number.\n\n"
  )

} else {

  stop(
    "Could not match taxonomy IDs to DADA2 ASVs.\n",
    "Expected IDs such as ASV1, ASV2, ASV3..."
  )
}


# ------------------------------------------------------------
# 7. Assign genus to every ASV
# ------------------------------------------------------------

genus_for_asv <- rep(
  "Unclassified",
  length(asv_names)
)

valid <- (
  !is.na(match_index) &
  match_index >= 1 &
  match_index <= length(genus)
)

genus_for_asv[valid] <-
  genus[match_index[valid]]


# Remove optional g__ prefix (before cleaning, so "g__unidentified" is caught)

genus_for_asv <- trimws(
  genus_for_asv
)

genus_for_asv <- sub(
  "^g__",
  "",
  genus_for_asv,
  ignore.case = TRUE
)


# Clean genus names

genus_for_asv[
  is.na(genus_for_asv) |
  genus_for_asv == "" |
  grepl("^unidentified", genus_for_asv, ignore.case = TRUE) |
  tolower(genus_for_asv) %in% c(
    "na",
    "nan",
    "none",
    "unknown",
    "unclassified",
    "uncultured",
    "no_match",
    "incertae_sedis",
    "incertae sedis"
  )
] <- "Unclassified"

cat(
  "ASVs with genus assignment: ",
  sum(genus_for_asv != "Unclassified"),
  " / ",
  length(genus_for_asv),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Combine ASVs into genera
# ------------------------------------------------------------

genus_levels <- unique(
  genus_for_asv
)

genus_table <- matrix(
  0,
  nrow = nrow(seqtab),
  ncol = length(genus_levels)
)

rownames(genus_table) <- rownames(seqtab)
colnames(genus_table) <- genus_levels


for (i in seq_along(genus_levels)) {

  selected <- which(
    genus_for_asv == genus_levels[i]
  )

  if (length(selected) == 1) {

    genus_table[, i] <-
      seqtab[, selected]

  } else {

    genus_table[, i] <-
      rowSums(
        seqtab[, selected, drop = FALSE]
      )
  }
}


genus_table <- as.data.frame(
  genus_table,
  check.names = FALSE
)

if (REMOVE_UNCLASSIFIED) {
  drop_cols <- tolower(names(genus_table)) %in% c("unidentified", "unclassified")
  if (any(drop_cols)) {
    cat("Removing unidentified/unclassified:",
        sum(genus_table[, drop_cols, drop = FALSE]), "reads\n\n")
    genus_table <- genus_table[, !drop_cols, drop = FALSE]
  }
  if (ncol(genus_table) == 0) stop("No classified genera left.")
}
# Save genus abundance table

write.table(
  genus_table,
  file = file.path(
    RESULTS_DIR,
    "genus_abundance_table.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)


# ============================================================
# 9. PLOT YOU REQUESTED
# ============================================================
#
# Total abundance of each fungal genus across all samples.
#
# The genera are ordered from most abundant to least abundant.
# ============================================================

total_abundance <- colSums(
  genus_table,
  na.rm = TRUE
)

genus_total <- data.frame(
  genus = names(total_abundance),
  Abundance = as.numeric(
    total_abundance
  ),
  stringsAsFactors = FALSE
)


# Sort from largest to smallest

genus_total <- genus_total[
  order(
    -genus_total$Abundance
  ),
]


# Show the top 25 genera

TOP_N <- 25

plot_total <- head(
  genus_total,
  min(
    TOP_N,
    nrow(genus_total)
  )
)


# Preserve abundance order

plot_total$genus <- factor(
  plot_total$genus,
  levels = plot_total$genus
)


# Create plot

p_total <- ggplot(
  plot_total,
  aes(
    x = genus,
    y = Abundance
  )
) +

  geom_col(
    width = 0.75
  ) +

  labs(
    title = "Fungal genus abundance",

    subtitle = paste0(
      "Top ",
      nrow(plot_total),
      " genera across all samples"
    ),

    x = "genus",

    y = "Abundance"
  ) +

  theme_minimal(
    base_size = 12
  ) +

  theme(

    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),

    panel.grid.major.x =
      element_blank(),

    panel.grid.minor =
      element_blank(),

    plot.title =
      element_text(
        face = "bold"
      )
  )


# Save plot

ggsave(
  filename = file.path(
    PLOTS_DIR,
    "genus_total_abundance.png"
  ),

  plot = p_total,

  width = 12,

  height = 7,

  dpi = 300
)


# ============================================================
# 10. RELATIVE ABUNDANCE PLOT
# ============================================================

sample_totals <- rowSums(
  genus_table
)

relative_table <- genus_table


for (i in seq_len(
  nrow(relative_table)
)) {

  if (sample_totals[i] > 0) {

    relative_table[i, ] <-
      relative_table[i, ] /
      sample_totals[i]
  }
}


# Convert to long format

relative_long <- do.call(
  rbind,

  lapply(
    seq_len(
      nrow(relative_table)
    ),

    function(i) {

      data.frame(

        Sample =
          rownames(relative_table)[i],

        genus =
          names(relative_table),

        RelativeAbundance =
          as.numeric(
            relative_table[i, ]
          ),

        stringsAsFactors = FALSE
      )
    }
  )
)


# Keep top genera separate

top_genera <- head(
  genus_total$genus,
  min(
    TOP_N,
    nrow(genus_total)
  )
)


relative_long$genus <-
  ifelse(
    relative_long$genus %in%
      top_genera,

    relative_long$genus,

    "Other"
  )


# Combine Other

relative_long <- aggregate(
  RelativeAbundance ~ Sample + genus,
  data = relative_long,
  FUN = sum
)


# Plot

p_relative <- ggplot(
  relative_long,
  aes(
    x = Sample,
    y = RelativeAbundance,
    fill = genus
  )
) +

  geom_col(
    width = 0.9
  ) +

  scale_y_continuous(
    labels = function(x) {
      paste0(
        round(x * 100),
        "%"
      )
    }
  ) +

  labs(

    title =
      "Fungal genus relative abundance",

    subtitle =
      "Top genera shown separately",

    x = "Sample",

    y = "Relative abundance",

    fill = "genus"
  ) +

  theme_minimal(
    base_size = 11
  ) +

  theme(

    axis.text.x =
      element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),

    panel.grid.major.x =
      element_blank(),

    panel.grid.minor =
      element_blank(),

    plot.title =
      element_text(
        face = "bold"
      )
  )


# Save

ggsave(
  filename = file.path(
    PLOTS_DIR,
    "genus_relative_abundance.png"
  ),

  plot = p_relative,

  width = 14,

  height = 8,

  dpi = 300
)


# ------------------------------------------------------------
# 11. BRAY-CURTIS PCoA + CLUSTERING
# ------------------------------------------------------------

sample_totals_for_pcoa <- rowSums(
  genus_table,
  na.rm = TRUE
)

keep_samples <- sample_totals_for_pcoa > 0

genus_table_pcoa <- genus_table[
  keep_samples,
  ,
  drop = FALSE
]

cat(
  "Samples used for Bray-Curtis PCoA: ",
  nrow(genus_table_pcoa),
  " / ",
  nrow(genus_table),
  "\n",
  sep = ""
)

if (nrow(genus_table_pcoa) < 3) {
  stop("Need at least 3 samples for clustering.")
}

# Bray-Curtis distance
bray <- vegdist(
  genus_table_pcoa,
  method = "bray"
)

# PCoA
pcoa <- cmdscale(
  bray,
  k = 2,
  eig = TRUE,
  add = TRUE
)

pcoa_points <- data.frame(
  Sample = rownames(genus_table_pcoa),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# CLUSTER SAMPLES
# ------------------------------------------------------------

# Change this number if you want more/fewer clusters
N_CLUSTERS <- 3

set.seed(123)

cluster_result <- kmeans(
  pcoa_points[, c("PCoA1", "PCoA2")],
  centers = N_CLUSTERS
)

pcoa_points$Cluster <- factor(
  cluster_result$cluster
)

# ------------------------------------------------------------
# PERCENTAGE OF VARIATION
# ------------------------------------------------------------

positive_eig <- pcoa$eig[
  pcoa$eig > 0
]

if (length(positive_eig) >= 2) {

  percent1 <- round(
    100 * positive_eig[1] / sum(positive_eig),
    1
  )

  percent2 <- round(
    100 * positive_eig[2] / sum(positive_eig),
    1
  )

} else {

  percent1 <- NA
  percent2 <- NA
}

# ------------------------------------------------------------
# PLOT
# ------------------------------------------------------------

p_pcoa <- ggplot(
  pcoa_points,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = Cluster
  )
) +

  geom_point(
    size = 3
  ) +

  labs(
    title = "Bray-Curtis PCoA",
    subtitle = paste0(
      "Samples clustered into ",
      N_CLUSTERS,
      " groups"
    ),
    x = paste0(
      "PCoA1",
      ifelse(
        is.na(percent1),
        "",
        paste0(" (", percent1, "%)")
      )
    ),
    y = paste0(
      "PCoA2",
      ifelse(
        is.na(percent2),
        "",
        paste0(" (", percent2, "%)")
      )
    ),
    color = "Cluster"
  ) +

  theme_minimal(
    base_size = 12
  ) +

  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# Save
ggsave(
  filename = file.path(
    PLOTS_DIR,
    "bray_curtis_pcoa_clustered.png"
  ),
  plot = p_pcoa,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# 12. Export taxonomy table
# ============================================================

taxonomy_export <- data.frame(

  ASV = asv_names,

  genus = genus_for_asv,

  stringsAsFactors = FALSE
)


write.table(

  taxonomy_export,

  file = file.path(
    RESULTS_DIR,
    "taxonomy_table.tsv"
  ),

  sep = "\t",

  quote = FALSE,

  row.names = FALSE
)


# ============================================================
# 13. Finished
# ============================================================

cat("\n")
cat("========================================\n")
cat("VISUALIZATION COMPLETE\n")
cat("========================================\n\n")

cat("Created:\n")

cat(
  "  results/plots/genus_total_abundance.png\n"
)

cat(
  "  results/plots/genus_relative_abundance.png\n"
)

cat(
  "  results/plots/bray_curtis_pcoa.png\n"
)

cat(
  "  results/genus_abundance_table.tsv\n"
)

cat(
  "  results/taxonomy_table.tsv\n"
)

cat("\nTop genera:\n")

print(
  head(
    genus_total,
    20
  )
)