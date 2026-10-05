# 03 ITS2 Taxonomy with dnabarcoder

In this demo we will be using dnabarcoder to do the taxonomy of our ASV table.

The reference sequences come from the UNITE ITS reference files

From the project root, you can choose which reference data from UNITE `dnabarcoder` will use.



```powershell
# Default (ITS2)
python scripts/03_taxonomy.py

# Select ITS1 region
python scripts/03_taxonomy.py --region its1

# Select full ITS region with custom CPU threads
python scripts/03_taxonomy.py -r its --ncpus 8
```

