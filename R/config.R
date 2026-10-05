# ============================================================
# R/config.R
#
# *** PLACEHOLDER — NOT YOUR ORIGINAL FILE ***
#
# Your original R/config.R was not among the files you gave me, so I could
# not include the real one without guessing at its contents. I looked at
# what scripts/02_dada2.R actually reads from this file, and the only
# value it uses is MIN_OVERLAP (passed to dada2::mergePairs()).
#
# 12 is dada2's own default for minOverlap, so this placeholder will let
# the pipeline run, but if your real config.R set this (or anything else)
# to a different value, replace this file with your original before
# trusting the results.
# ============================================================

MIN_OVERLAP <- 12
