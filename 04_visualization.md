# 04 Visualization

This step reads:

```text
results/seqtab_nochim.rds
results/asvs.merged.classification
```

and creates the main R outputs described by the supplied tutorial:
* taxonomic composition bar plot
* Bray–Curtis PCoA
* exported combined ASV/taxonomy table

You should define your reference data, and make sure it is the same as the previous step. You can also choose which taxonomy you will visualize

Example usage:

```powershell
# Default: ITS2 region at the Genus level
Rscript scripts/04_visualization.R

# Specify region (its1, its2, or its)
Rscript scripts/04_visualization.R --region its1

# Specify taxonomic rank (phylum, class, order, family, genus)
Rscript scripts/04_visualization.R --rank phylum

# Combine both arguments
Rscript scripts/04_visualization.R --region its1 --rank family

# Shorthand flags
Rscript scripts/04_visualization.R -r its -k class
```