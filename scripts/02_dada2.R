# ============================================================
# 02_dada2.R
# Paired-end fungal ITS2 DADA2 pipeline for Windows
#
# Based on the supplied DADA2 tutorial.
# Corrupted or incomplete paired FASTQ samples are skipped.
#
# Expected input:
#   data/fastq/SAMPLE_1.fastq.gz
#   data/fastq/SAMPLE_2.fastq.gz
#
# Run from the project root:
#   Rscript scripts/02_dada2.R
# ============================================================

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

source(file.path(PROJECT_ROOT, "R", "config.R"))

suppressPackageStartupMessages(library(dada2))

FASTQ_DIR <- file.path(PROJECT_ROOT, "data", "fastq")
FILTERED_DIR <- file.path(PROJECT_ROOT, "data", "filtered")
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
DNABARCODER_DIR <- file.path(PROJECT_ROOT, "dnabarcoder_input")

dir.create(FILTERED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DNABARCODER_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Find paired FASTQ files
# ------------------------------------------------------------

fnFs.all <- sort(list.files(
  FASTQ_DIR,
  pattern = "_1\\.fastq\\.gz$",
  full.names = TRUE
))

fnRs.all <- sort(list.files(
  FASTQ_DIR,
  pattern = "_2\\.fastq\\.gz$",
  full.names = TRUE
))

sample.names.F <- sub("_1\\.fastq\\.gz$", "", basename(fnFs.all))
sample.names.R <- sub("_2\\.fastq\\.gz$", "", basename(fnRs.all))

common.samples <- intersect(sample.names.F, sample.names.R)

if (length(common.samples) == 0) {
  stop(
    "No matching paired FASTQ files were found in: ",
    FASTQ_DIR
  )
}

fnFs <- file.path(
  FASTQ_DIR,
  paste0(common.samples, "_1.fastq.gz")
)

fnRs <- file.path(
  FASTQ_DIR,
  paste0(common.samples, "_2.fastq.gz")
)

names(fnFs) <- common.samples
names(fnRs) <- common.samples

cat("\nFound ", length(common.samples), " candidate paired samples.\n", sep = "")

# ------------------------------------------------------------
# 2. Check gzip/FASTQ integrity
# ------------------------------------------------------------
#
# A sample is skipped if either mate is incomplete, corrupted,
# unreadable, or contains invalid FASTQ records.
# ------------------------------------------------------------

check_fastq <- function(path) {
  result <- tryCatch({
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)

    reads <- 0L

    repeat {
      lines <- readLines(con, n = 4L)

      if (length(lines) == 0L) {
        break
      }

      if (length(lines) != 4L) {
        stop("Incomplete FASTQ record")
      }

      if (
        !startsWith(lines[1], "@") ||
        !startsWith(lines[3], "+")
      ) {
        stop("Invalid FASTQ record")
      }

      if (nchar(lines[2]) != nchar(lines[4])) {
        stop("Sequence and quality lengths differ")
      }

      reads <- reads + 1L
    }

    list(ok = reads > 0L, reads = reads, error = NA_character_)
  }, error = function(e) {
    list(ok = FALSE, reads = 0L, error = conditionMessage(e))
  })

  result
}

integrity <- data.frame(
  sample = common.samples,
  forward_ok = FALSE,
  reverse_ok = FALSE,
  forward_reads = 0L,
  reverse_reads = 0L,
  status = "SKIP",
  reason = "",
  stringsAsFactors = FALSE
)

cat("\nChecking FASTQ integrity...\n")

for (i in seq_along(common.samples)) {
  sample <- common.samples[i]

  f <- check_fastq(fnFs[i])
  r <- check_fastq(fnRs[i])

  integrity$forward_ok[i] <- f$ok
  integrity$reverse_ok[i] <- r$ok
  integrity$forward_reads[i] <- f$reads
  integrity$reverse_reads[i] <- r$reads

  if (f$ok && r$ok) {
    integrity$status[i] <- "KEEP"
    integrity$reason[i] <- "Both mates valid"
    cat("KEEP  ", sample, "\n", sep = "")
  } else {
    reasons <- c()

    if (!f$ok) {
      reasons <- c(reasons, paste0("forward: ", f$error))
    }

    if (!r$ok) {
      reasons <- c(reasons, paste0("reverse: ", r$error))
    }

    integrity$reason[i] <- paste(reasons, collapse = " | ")

    cat(
      "SKIP  ", sample,
      " - ", integrity$reason[i],
      "\n",
      sep = ""
    )
  }
}

write.table(
  integrity,
  file = file.path(RESULTS_DIR, "fastq_integrity_report.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

keep <- integrity$status == "KEEP"

if (!any(keep)) {
  stop(
    "No complete paired FASTQ samples remain. ",
    "See results/fastq_integrity_report.tsv."
  )
}

sample.names <- common.samples[keep]
fnFs <- fnFs[keep]
fnRs <- fnRs[keep]

# ------------------------------------------------------------
# 3. Define filtered FASTQ paths
# ------------------------------------------------------------

filtFs <- file.path(
  FILTERED_DIR,
  paste0(sample.names, "_1_filtered.fastq.gz")
)

filtRs <- file.path(
  FILTERED_DIR,
  paste0(sample.names, "_2_filtered.fastq.gz")
)

names(filtFs) <- sample.names
names(filtRs) <- sample.names

# ------------------------------------------------------------
# 4. Filter and trim
# ------------------------------------------------------------
#
# For ITS, no fixed truncLen is used because ITS lengths vary.
# Windows: multithread = FALSE.
# ------------------------------------------------------------

cat("\nFiltering ", length(sample.names), " complete pairs...\n", sep = "")

out <- filterAndTrim(
  fnFs,
  filtFs,
  fnRs,
  filtRs,
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE
)

print(out)

keep.filtered <- out[, "reads.out"] > 0

if (!any(keep.filtered)) {
  stop("No reads passed filtering for any valid paired sample.")
}

if (any(!keep.filtered)) {
  cat(
    "\nSamples with zero reads after filtering:\n",
    paste(sample.names[!keep.filtered], collapse = "\n"),
    "\n"
  )
}

sample.names <- sample.names[keep.filtered]
fnFs <- fnFs[keep.filtered]
fnRs <- fnRs[keep.filtered]
filtFs <- filtFs[keep.filtered]
filtRs <- filtRs[keep.filtered]
out <- out[keep.filtered, , drop = FALSE]

# ------------------------------------------------------------
# 5. Learn error rates
# ------------------------------------------------------------

cat("\nLearning forward error rates without quality scores...\n")
errF <- learnErrors(
  filtFs,
  errorEstimationFunction = noqualErrfun,
  multithread = FALSE
)

cat("\nLearning reverse error rates without quality scores...\n")
errR <- learnErrors(
  filtRs,
  errorEstimationFunction = noqualErrfun,
  multithread = FALSE
)

# ------------------------------------------------------------
# 6. Denoise reads
# ------------------------------------------------------------

cat("\nDenoising forward reads...\n")
dadaFs <- dada(filtFs, err = errF, multithread = FALSE)

cat("\nDenoising reverse reads...\n")
dadaRs <- dada(filtRs, err = errR, multithread = FALSE)

# ------------------------------------------------------------
# 7. Merge paired reads
# ------------------------------------------------------------

cat("\nMerging paired reads...\n")

mergers <- mergePairs(
  dadaFs,
  filtFs,
  dadaRs,
  filtRs,
  minOverlap = MIN_OVERLAP,
  verbose = TRUE
)

names(mergers) <- sample.names

# ------------------------------------------------------------
# 8. Construct sequence table
# ------------------------------------------------------------

seqtab <- makeSequenceTable(mergers)

cat(
  "\nSequence table: ",
  nrow(seqtab),
  " samples x ",
  ncol(seqtab),
  " ASVs\n",
  sep = ""
)

cat("\nSequence length distribution:\n")
print(table(nchar(getSequences(seqtab))))

# ------------------------------------------------------------
# 9. Remove chimeras
# ------------------------------------------------------------

cat("\nRemoving chimeras...\n")

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = FALSE,
  verbose = TRUE
)

# ------------------------------------------------------------
# 10. Track reads
# ------------------------------------------------------------

getN <- function(x) {
  sum(getUniques(x))
}

track <- cbind(
  out,
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)

colnames(track) <- c(
  "input",
  "filtered",
  "denoisedF",
  "denoisedR",
  "merged",
  "nonchim"
)

rownames(track) <- sample.names

write.table(
  track,
  file = file.path(RESULTS_DIR, "read_tracking.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

print(track)

# ------------------------------------------------------------
# 11. Save outputs
# ------------------------------------------------------------

saveRDS(
  seqtab.nochim,
  file = file.path(RESULTS_DIR, "seqtab_nochim.rds")
)

write.table(
  seqtab.nochim,
  file = file.path(RESULTS_DIR, "asv_table.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

asv.seqs <- getSequences(seqtab.nochim)
asv.fasta <- file.path(DNABARCODER_DIR, "asvs.fasta")

fasta.con <- file(asv.fasta, open = "wt")

for (i in seq_along(asv.seqs)) {
  writeLines(paste0(">ASV", i), fasta.con)
  writeLines(asv.seqs[i], fasta.con)
}

close(fasta.con)

# ------------------------------------------------------------
# 12. Final summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("DADA2 PIPELINE COMPLETE\n")
cat("============================================\n")

cat("Samples processed: ", nrow(seqtab.nochim), "\n", sep = "")
cat("ASVs after chimera removal: ", ncol(seqtab.nochim), "\n", sep = "")

cat("\nReports and outputs:\n")
cat("  ", file.path(RESULTS_DIR, "fastq_integrity_report.tsv"), "\n", sep = "")
cat("  ", file.path(RESULTS_DIR, "read_tracking.tsv"), "\n", sep = "")
cat("  ", file.path(RESULTS_DIR, "seqtab_nochim.rds"), "\n", sep = "")
cat("  ", file.path(RESULTS_DIR, "asv_table.tsv"), "\n", sep = "")
cat("  ", asv.fasta, "\n", sep = "")
