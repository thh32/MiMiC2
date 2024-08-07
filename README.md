<img src="https://github.com/thh32/MiMiC2/assets/" alt="MiMiC2-logo" width="555" height="128"/>

## Installation Instructions
1. Clone the repository.  
```bash
git clone https://github.com/Matt-Schmitz/InvestiGut.git
```
2. Enter the InvestiGut folder.  
```bash
cd InvestiGut
```
3. Create the environment with [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).  
```bash
mamba create --no-channel-priority -n investigut \
    -c bioconda -c conda-forge \
    "python=3.11" "numpy=1.24.3" "scipy=1.10.1" \
    "conda-forge::matplotlib-base" "seaborn=0.13.0" \
    "pandas=1.5.3" "statsmodels=0.13.5" "ete3=3.1.2" \
    "openpyxl=3.0.10" "bioconda::diamond=2.1.8"
```
4. Activate the environment.  
```bash
mamba activate investigut
```
5. Download proteins and create DIAMOND database.  After the DIAMOND database has been created, the compressed fasta data can be deleted.

> [!WARNING]
> The compressed data take up 79GB (two files).
> The DIAMOND database will take up 311GB.

```bash
./download_fasta.sh
```
6. Add the InvestiGut folder to your `~/.bashrc` and apply the changes.  
- Add to .bashrc: `export PATH="/PATH/TO/InvestiGut:$PATH"`
- Enter in terminal:
```bash
source ~/.bashrc
```
- Reactivate the mamba environment:
```bash
mamba activate investigut
```

7. Run InvestiGut on chosen proteins.  
Basic usage:

```bash
investigut.py -i /PATH/TO/FASTA/FILE
```













pandas

seaborn

tqdm


glpk


bacarena

sudo apt-get -y install libglpk-dev


R
Rscript -e 'remotes::install_github("SysBioChalmers/sybil")'
Rscript -e 'remotes::install_github("euba/bacarena")'












  
Options list:

Any arguments not recognized by InvesitGut will be passed on to DIAMOND (e.g., `--threads`, `--faster`, `--ultra-sensitive`).
```
-h, --help            show this help message and exit
-i {INPUT}            input fasta file
-o {OUTPUT}           output folder
-s                    single mode (default)
-m                    multi mode
-d {DMND_TSV_OUTPUT}  optional DIAMOND .tsv output file to skip alignment
--query-cover {NUM}   DIAMOND query-coverage (default 90.0)
--subject-cover {NUM}
                      DIAMOND subject-coverage (default 90.0)
--id {NUM}            DIAMOND %id (default 90)
--low                 overrides default DIAMOND query/sujbect coverage and %id giving all a
                      value of 50
```
## Example Usage

In the example folder, a file called "seaweed.fa" contains two bacterial proteins, Bp1670 and Bp1689, that are involved in seaweed digestion (taken from https://www.nature.com/articles/nature08937). In multi mode (`-m`), metagenomes and metagenome-assembled-genomes (MAGs) containing all proteins in the fasta file are counted as positive, whereas in single mode (`-s`), metagenomes and MAGs are separately examined for each protein sequence. The following command searches for matches to both proteins (`-m`) using a subject/query-coverage of 50% and a percetage identity threshhold of 50% (`--low`). 
```bash
investigut.py -i /PATH/TO/seaweed.fa -m --low
```
The resulting output can be found under `./examples/Bp1670 + Bp1689`. In this folder, the "overview.txt" file contains global prevalence within the metagenomes, disease prevalence statistics, country prevalence statistics, and other demographic factor data (smokers vs. non-smokers, BMI, gender, age by decade of life, and antibiotic usage). The MAG overview data include a list of positive gut bacteria, occurrence of the protein(s) by taxonomic rank, and the cumulative relative abundance of species containing the protein(s) of interest within the two cohorts of the origin data (https://www.nature.com/articles/s41467-022-31502-1). 

The following is an example of an auto-generated figure of prevalence by country found within the output folder after running the seaweed digestion proteins mentioned above through InvestiGut.

![Prevalence_of_Bp1670 + Bp1689_by_Country](https://github.com/Matt-Schmitz/InvestiGut/assets/34464190/12cb4a4a-4cd5-47cd-af14-803e949c310c)
