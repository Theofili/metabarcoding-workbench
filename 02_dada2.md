# 02 DADA2 Processing

In this step we use **DADA2** to turn raw paired-end sequencing reads into a table of **ASVs** (amplicon sequence variants): exact, error-corrected sequences that represent the distinct organisms in your samples. The ASV sequences are exported as a FASTA file for the taxonomy step.

## What the script does

1.  **Quality filtering.** Reads with low-quality bases, too many expected errors, or ambiguous bases (N) are trimmed or discarded, so poor data doesn't create false variants.

2.  **Error learning.** DADA2 learns the error pattern of your sequencing run (how often each base is miscalled as another). Forward and reverse reads are modeled separately.

3.  **DADA2 denoising.** Using the learned error model, DADA2 decides whether each rare sequence is a real biological variant or just a sequencing error of a more abundant one, and corrects the errors. This is what produces exact ASVs rather than fuzzy OTU clusters.

4.  **Paired-read merging.** Each forward read is merged with its matching reverse read where they overlap, giving the full amplicon sequence. Pairs that don't overlap cleanly are dropped.

5.  **Chimera removal.** Chimeras (artificial sequences formed when two different templates fuse during PCR) are detected and removed.

6.  **ASV FASTA export.** The final ASV sequences are written to `dnabarcoder_input/asvs.fasta`, which is the input for step 03.

## Run

```         
Rscript scripts\02_dada2.R
```

**This script might take a while to finish**, especially the error-learning and denoising steps, which scale with the number of samples and reads.