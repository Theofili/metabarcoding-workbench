# 03 ITS Taxonomy with dnabarcoder

In this step we assign a taxonomic identity (kingdom → species) to each ASV produced by the DADA2 step, using **dnabarcoder** and the **UNITE** ITS reference database.

## What the script does

1.  **Checks inputs.** Confirms that `asvs.fasta`, the UNITE reference files for your chosen region, dnabarcoder, and BLAST+ (`blastn`, `makeblastdb`) are all available. If anything is missing you get a clear message before any slow work starts.

2.  **BLAST best-match search.** Each ASV is BLASTed against the UNITE reference sequences, and only the single best hit per ASV is kept. Very short alignments (under 50 bp) are penalized so a tiny perfect match can't beat a longer, slightly less identical one. A BLAST database is built on the first run and reused afterwards.

3.  **dnabarcoder classification.** The best hits are passed to dnabarcoder, which uses *similarity cutoffs for each taxonomic rank* to decide how deep the assignment can be trusted. An ASV with a very close match may be named to species; one with a distant match may only be assigned to genus, family, or higher.

## Reference files required

For each region, these three files must be in the `references/` folder:

| File | Purpose |
|----|----|
| `unite2025<REGION>.fasta` | Reference sequences (used to build the BLAST database) |
| `unite2025<REGION>.classification` | Taxonomic lineage of each reference sequence |
| `unite2025<REGION>.unique.cutoffs.best.json` | Similarity threshold for each taxonomic rank |

`<REGION>` is `ITS1`, `ITS2`, or `ITS`.

## Run

From the project root, choose which UNITE reference region to use:

```         
# Default (ITS2) 
python scripts/03_taxonomy.py  

# Select ITS1 region 
python scripts/03_taxonomy.py --region its1 

# Select full ITS region with custom CPU threads 
python scripts/03_taxonomy.py -r its --ncpus 8
```

| Option | Meaning |
|----|----|
| `--region`, `-r` | `its1`, `its2`, or `its` (default `its2`). Pick the region that matches the primers you amplified. |
| `--ncpus` | Number of CPU threads for BLAST (default 4). |

## Output

Results are written to `results/dnabarcoder/`:

- `asvs.<reference>_BLAST.bestmatch`: best BLAST hit for every ASV (reference ID, score, similarity, coverage)

- `asvs.<reference>_BLAST.classified`: **final taxonomy for each ASV**