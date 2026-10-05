# 01 Fetch and Prepare the Dataset

The demo tutorial uses BioProject:

**PRJNA356769**

## Paper Details

| Vu D, de Vries M, Gerrits van den Ende B, et al. (2026). “Advancing Yeast Identification Using High-Throughput DNA Barcode Data From a Curated Culture Collection.” Molecular Ecology Resources 26(1): e70082.https://doi.org/10.1111/1755-0998.70082

*That paper’s metabarcoding demonstration reclassifies the Human Microbiome Project gut mycobiome dataset from:*

| Nash AK, Auchtung TA, Wong MC, et al. (2017) “The gut mycobiome of the Human Microbiome Project healthy cohort.” Microbiome 5:153. https://doi.org/10.1186/s40168-017-0373-4

## Data Details

- Primary Target (ITS2): Amplified and sequence the Internal Transcribed Spacer 2 region of the eukaryotic rRNA operon. The universal fungal barcode.

- Sample count: Analyzed 317 stool samples collected longitudinally from healthy adult volunteers in the Human Microbiome Project (HMP) cohort.

- Repository Accession: the raw sequence data is deposited in the NCBI Sequence Read Archive (SRA) under BioProject: PRJNA356769.

- The dataset is heavily dominated by yeast species. Eight of the most abundant genera are yeast, led by:

  - Saccharomyces (96.8% of samples)

  - Malassezia (88.3% of samples)

  - Candida (80.8% of the samples)

Run:

```powershell
python scripts\01_download.py PRJNA356769
```

If you want to install some of the samples from the study use:

```powershell
python scripts\01_download.py PRJNA356769 --limit 5 # wwill donwload just the first five sample
```
The project contains multiple marker types, including 18S, so the important part is to select the ITS2 paired-end reads rather than downloading every FASTQ file.
 When running this command, you access NCBI's SRA database and start downloading fastq.gz files related to the BioProject Number. If you change the command line argument BioProject number you will donwload related data

*This step can be skipped if you have your own data uploaded to the corresponding folder `metabarcoding-workbench/data/fastq/`.*