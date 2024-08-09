<img src="https://github.com/thh32/MiMiC2/images/2024-08-08_Logo.png" alt="MiMiC2-logo" width="555" height="128"/>

## Installation Instructions
1. Clone the repository.  
```bash
git clone https://github.com/thh32/MiMiC2.git
```

2. Enter the MiMiC2 folder.  
```bash
cd MiMiC2
```

3. Create the environment with [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).  
```bash
mamba create --no-channel-priority -n mimic2 \
    -c bioconda -c conda-forge \
    "python=3.11" "numpy=1.24.3" "scipy=1.10.1" \
    "conda-forge::matplotlib-base" "seaborn=0.13.0" \
    "pandas=1.5.3" "statsmodels=0.13.5" "tdqm" \
    "openpyxl=3.0.10" "bioconda::diamond=2.1.8"
```

4. Activate the environment.  
```bash
mamba activate mimic2
```

5. Install glpk and the required R packages.

> [!WARNING]
>  **Sudo access is required for this**
> If you do not have sudo access please talk to your system administrator.

Install glpk-dev:
```
sudo apt-get -y install libglpk-dev
```

Next, install the R packages required for metabolic modelling:
```
Rscript -e 'remotes::install_github("SysBioChalmers/sybil")'
```
```
Rscript -e 'remotes::install_github("euba/bacarena")'
```

6. Add the MiMiC2 folder to your `~/.bashrc` and apply the changes.  
- Add to .bashrc: `export PATH="/PATH/TO/MiMiC2:$PATH"`
- Enter in terminal:
```bash
source ~/.bashrc
```
- Reactivate the mamba environment:
```bash
mamba activate mimic2
```

7. Run MiMiC2 on your chosen genome collection and metagenomic samples:  
Basic usage:

```bash
MiMiC2.py -g /PATH/TO/GENOME-COLLECTION -t /PATH/TO/TAXONOMIC-FILE --taxonomiclevel s -s /PATH/TO/SAMPLES -m /PATH/TO/METADATA --group GROUP --models /PATH/TO/MODELS/FOLDER -c 10 -o OUTPUT-PREFIX
```

  
Options list:
```
  -h, --help            show this help message and exit
  -s {INPUT}, --samples {INPUT}
                        Pfam vector file of all metagenomic samples to be
                        studied.
  -m {INPUT}, --metadata {INPUT}
                        Metadata file detailing the group assignment of each
                        sample.
  --group {INPUT}       Name of the group of interest for SynCom creation.
  -p {INPUT}, --pfam {INPUT}
                        Pfam file e.g. XXXXXXXX.
  -g {INPUT}, --genomes {INPUT}
                        Pfam vector file of genome collection.
  --models {INPUT}      Folder containing metabolic models for each genome.
                        Must be provided as RDS files such as those provided
                        by GapSeq.
  -t {INPUT}, --taxonomy {INPUT}
                        Taxonomic assignment of each genome.
  --taxonomiclevel {INPUT}
                        Taxonomic level for filtering (species = s,genus = g,
                        class = c, order = o, phyla = p).
  -c {INT}, --consortiasize {INT}
                        Define the SynCom size the user is after.
  --corebias {FLOAT}    The additional weighting provided to Pfams core to the
                        studied group (default = 0.0005).
  --groupbias {FLOAT}   The additional weighting provided to Pfams
                        significantly enriched in the studied group (default =
                        0.0012).
  --prevfilt {FLOAT}    Prevalence filtering threshold for shortlisting
                        genomes for inclusion in the final SynCom selection
                        (default = 33.3).
  -o {OUTPUT}, --output {OUTPUT}
                        Prefix for all output files e.g. HuSynCom.
  --iterations {INT}    Change the number of iterations to select sample
                        specific strains in step 1.
  --exclusion {INPUT}   Provide file which includes a csv list of genome names
                        to be excluded during SynCom selection.
```


# Datasets Provided
## Genome Collection Datasets
In the `datasets/isolate_collections/` folder we already provide preprocessed input files for a publicly available isolate collections from mice (miBCII), ruminants (Hungate1000), pigs (PiBAC), and the human gut (HiBC).  

We have also processed two large collections which include MAGs of uncultured taxa in `datasets/mag_collections/`, preventing experimental use of these SynComs. These are provided to allow initial study of SynComs prior to isolation, allowing researchers to target the isolation towards microbes of particular interest.

For each of these, the Pfam profiles and taxonomic assignments files are provided, allowing the direct use without further processing.

## Environmental Datasets
In the `Datasets/Environmental` folder we provide Pfam profiles for all samples studied within the MiMiC2 paper. 

## Example Run of MiMiC2
As an example of how MiMiC2 can be used, we will repeat the creation of the IBD SynCom from the MiMiC2 paper.

Run the code below, which uses the HiBC collection, along with the IBDMDB sample collection, with GEMs made with GapSeq:
```bash
MiMiC2.py -g datasets/isolate_collections/HiBC/HiBC_profile.txt -t /PATH/TO/TAXONOMIC-FILE --taxonomiclevel s -s /PATH/TO/SAMPLES -m /PATH/TO/METADATA --group GROUP --models /PATH/TO/MODELS/FOLDER -c 10 -o Test-IBD-SynCom
```
The  output can be found under `./Test-IBD-SynCom/`. In this folder, the "overview.txt" file contains global prevalence within the metagenomes, disease prevalence statistics, country prevalence statistics, and other demographic factor data (smokers vs. non-smokers, BMI, gender, age by decade of life, and antibiotic usage). The MAG overview data include a list of positive gut bacteria, occurrence of the protein(s) by taxonomic rank, and the cumulative relative abundance of species containing the protein(s) of interest within the two cohorts of the origin data (https://www.nature.com/articles/s41467-022-31502-1). 
