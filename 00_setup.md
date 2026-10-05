# 00 Setup

For this workflow, you will need to have Python, R, BLAST+.

## 0. Navigate to the correct project folder

All files will be saved inside this folder.

**For this demo the working directory is:**
```text
D:\metabarcoding-workbench
```

## 1. Download the files..

The folder should contain:

```text
00_setup.md
01_download_data.md
02_dada2.md
03_taxonomy.md
04_visualization.md
requirments.txt
scripts/
references/
R/
```

## Check Python

Run:

```powershell
py --version
```

You should see:

```text
Python 3.14.7 #or newer
```

## Check R

Run in Powershell:
```powershell
$env:Path += ";C:\Program Files\R\R-4.4.1\bin"
R.exe --version
```
You should see:
```text
R version 4.4.1 
```
Then check Rscript:
```powershell
Rscript --version
```

You should see:

```text
Rscript (R) version 4.4.1
```

## 2. Create or use the Project enviroment

if `.venv` does not exist inside `metabarcoding-workbench`, create it:

```powershell
py -m venv .venv
```

You can activate it:

```powershell
.\.venv\Scripts\Activate.ps1 
```

## 3. Install Python Packages

```powershell
.\.venv\Scripts\python.exe -m pip install --upgrade pip 
```

Then:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt 
```

This installs:

* cutadapt
* Biopython
* scikit-learn
* matplotlib
* pandas

Check the Python setup:

```powershell
.\.venv\Scripts\python.exe scripts\00_check_setup.py 
```

*I don't have cutadapt for some reason now*

## 4. Install R packages

Use `Rscript` directly from the **PowerShell terminal**

Install `ggplot`, `vegan` and `BiocManager`:

```powershell
Rscript -e "install.packages(c('ggplot2','vegan','BiocManager'), repos='https://cloud.r-project.org')"  
```
Install DADA2 AND ShortRead:

```powershell
Rscript -e "BiocManager::install(c('dada2','ShortRead'), ask=FALSE, update=FALSE, force=TRUE)"
```
This step may take several minutes

## Check the R packages

```powershell
Rscript -e "library(ggplot2); library(vegan); library(dada2); library(ShortRead); cat('All required R packages loaded successfully.\n')" 
```
You want to see:
```text
All required R packages loaded successfully.
```

If any package has not been downloaded, **stop** and fix that package before continuing.


## 00a - Install BLAST+, dnabarcoder and UNITE reference files

For this demo dnabarcoder will be used to create the taxonomy files.The reference sequences come from the **UNITE+INSD 2025 fungal ITS2 dataset**, which was prepared for dnabarcoder.

## 1. Install NCBI BLAST+

From NCBI download the Windows installer:
```powershell
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```

Run to install:
```powershell
C:\Program Files\NCBI\blast-2.16.0+\bin
```
*The path may differ to different computers.*

Add that `bin` folder to PATH for the current terminal

```powershell
$env:Path += ";C:\Users\theos\Downloads\ncbi-blast-2.17.0+-x64-win64\ncbi-blast-2.17.0+\bin"
```
Verify:
```powershell
blastn -version
```

## 2. Get dnabarcoder

Using git:
```powershell
git clone https://github.com/vuthuyduong/dnabarcoder.git
```
You must have a new folder `metabarcoding-workbench/dnabarcoder`

## 3. Get the UNITE ITS2 reference files

Source: https://zenodo.org/records/22210122
(UNITE+INSD 2025 Fungal ITS, ITS1, and ITS2 Reference Sequences with their similarity cutoffs)

Download the zip (about 290 MB):

```powershell
Invoke-WebRequest -Uri "https://zenodo.org/records/22210122/files/UNITE_2025.zip?download=1" -OutFile UNITE_2025.zip
```

*This step may take a while*

Optional: check the download is complete. The MD5 should be `edc1af5987818dd0d3b2d9f664d918fa`:

```powershell
Get-FileHash UNITE_2025.zip -Algorithm MD5
```

Unzip it and list the ITS2 files:

*If your data come from **ITS** or **ITS1** you can also change the ITS2 argument and donwload the reference data suited for you.*


```powershell
Expand-Archive UNITE_2025.zip -DestinationPath unite_tmp
Get-ChildItem unite_tmp -Recurse -Filter "*ITS2*" | Select-Object Name, Length
```

You need three ITS2 files: the FASTA, the classification, and the cutoffs JSON.
For the 2024 release they were named:

```text
unite2024ITS2.fasta
unite2024ITS2.unique.classification
unite2024ITS2.unique.cutoffs.best.json
```

so expect the same pattern with `2025`. The `.unique` files are the smaller,
de-duplicated versions. **Use the names your listing actually shows.**

Copy them into `references\`:

```powershell
Copy-Item (Get-ChildItem unite_tmp -Recurse -Filter "unite2025ITS2.fasta").FullName                        references\
Copy-Item (Get-ChildItem unite_tmp -Recurse -Filter "unite2025ITS2.classification").FullName        references\
Copy-Item (Get-ChildItem unite_tmp -Recurse -Filter "unite2025ITS2.unique.cutoffs.best.json").FullName      references\
```

Check that all three arrived:

```powershell
Get-ChildItem references\ -Filter "unite*ITS2*"
```

Once copied you can delete the temporary files:

```powershell
Remove-Item unite_tmp -Recurse -Force
Remove-Item UNITE_2025.zip
```