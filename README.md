<img src="https://github.com/thh32/MiMiC2/blob/main/images/2024-08-08_Logo.png" alt="MiMiC2-logo" width="555" height="150"/>

## What is MiMiC2?

MiMiC2 is a bioinformatic pipeline for the selection of a few microbial genomes that functionally represent an entire ecosystem, termed a synthetic community (SynCom). 

If you want to install MiMiC2, go to the [intallation instructions](#installation-instructions).

If you want help running MiMiC2, go to the [running MiMiC2 guide](#running-mimic2).

If you want to understand all the options, go to the [options list](#options).

If you want premade datasets, look at our [datasets provided](#datasets-provided) section.

## MiMiC2 workflow
The general process can be seen in the workflow:

<p align="center">
<img src="https://github.com/thh32/MiMiC2/blob/main/images/2024-08-08_MiMiC2-workflow.png" alt="MiMiC2-workflow" width="350" height="750"/>
</p>


## Guide to run MiMiC2
MiMiC2 consists of a few major steps, but before that can begin you must prepare your data.

### Data Preparation

You must use the `MiMiC2-BUTLER.py` script to convert a folder containing multiple genomes/samples Pfam annotations into a single Pfam profile file of the entire dataset.
  
  `MiMiC2-BUTLER.py` needs a few bits of information to run, detailed under the [options list](#options) . We have provided example data to help you understand the process. The example of HiBC can be run using the code below:
  ```
  MiMiC2-BUTLER.py -s datasets/isolate_collections/HiBC/Pfams/ -p datasets/core/Pfam-A.clans.tsv -t hmmscan -e .hmmer -o HiBC-0.6-profile.txt
  ```


### Running MiMiC2
Once you have a Pfam profile of both your environment samples, and your genome collection, you can run MiMiC2. The options for running MiMiC2 are detailed in the [options list](#options) section. 
```bash
MiMiC2.py -g /PATH/TO/GENOME-COLLECTION -t /PATH/TO/TAXONOMIC-FILE --taxonomiclevel s -s /PATH/TO/SAMPLES -m /PATH/TO/METADATA --group GROUP --models /PATH/TO/MODELS/FOLDER -c 10 -o OUTPUT-PREFIX
```


### Example Run of MiMiC2
As an example of how MiMiC2 can be used, we will repeat the creation of the IBD SynCom from the MiMiC2 paper.

Run the code below, which uses the HiBC collection, along with the IBDMDB sample collection, with GEMs made with GapSeq:
```bash
MiMiC2.py -g datasets/isolate_collections/HiBC/HiBC_profile.txt -t /PATH/TO/TAXONOMIC-FILE --taxonomiclevel s -s /PATH/TO/SAMPLES -m /PATH/TO/METADATA --group GROUP --models /PATH/TO/MODELS/FOLDER -c 10 -o Test-IBD-SynCom
```
The  output can be found under `./Test-IBD-SynCom/`.


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
    "pandas=1.5.3" "tdqm" 
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



## Options
  
### MiMiC2-BUTLER.py options
```
  -h, --help            show this help message and exit
  -s {INPUT}, --samples {INPUT}
                        Provide a folder which contains all of your Pfam
                        annotated genomes/metagenomes.
  -p {INPUT}, --pfam {INPUT}
                        Pfam file e.g. Pfam-A.clans.csv, provided for Pfam v32
                        in `datasets/core/`
  -t {TEXT}, --tool {TEXT}
                        State the tool used to annotate the geomes against the
                        Pfam database: `hmmsearch` or `hmmscan`
  -o {OUTPUT}, --output {OUTPUT}
                        Prefix for all the Pfam-profile file e.g. HuSynCom.
  -e {TEXT}, --extension {TEXT}
                        Provide the extension for your Pfam annotation files.
```


### MiMiC2.py options
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
                        Pfam file e.g. Pfam-A.clans.csv, provided for Pfam v32
                        in `datasets/core/`
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
In the `datasets/isolate_collections/` folder we provide preprocessed input files for a publicly available isolate collections, allowing construction of any designed SynComs.  

We have also processed two large collections which include MAGs of uncultured taxa in `datasets/mag_collections/`, preventing experimental use of these SynComs. These are provided to allow initial study of SynComs prior to isolation, allowing researchers to target the isolation towards microbes of particular interest.

For each of these, the Pfam profiles and taxonomic assignments files are provided, allowing the direct use without further processing.

|Collection name | Environment | Isolates or Genomes|
|----------------|-----------|------------|
| [PiBAC](https://www.nature.com/articles/s41467-020-19929-w) | Pig gut | Isolates |
| [HiBC](https://www.biorxiv.org/content/10.1101/2024.06.20.599854v1) | Human gut | Isolates |
| [miBCII](https://www.sciencedirect.com/science/article/pii/S193131282200467X) | Mouse gut | Isolates |
| [Hungate1000](https://www.nature.com/articles/nbt.4110) | Rumen | Isolates |
| [GTDB r202](https://academic.oup.com/nar/article/50/D1/D785/6370255) | N/A | Genomes |
| [Pasolli et al 2016](https://www.sciencedirect.com/science/article/pii/S0092867419300017) | Human microbiome | Genomes |



## Environmental Datasets
In the `Datasets/Environmental` folder we provide Pfam profiles for all samples studied within the MiMiC2 paper. 


|Study | Environment | Condition|
|----------------|-----------|------------|
|[Shabat et al, 2016.](https://doi.org/10.1038/ismej.2016.62) | Bovine rumen | N/A|
| [Lesker et al, 2020](https://doi.org/10.1016/j.celrep.2020.02.036) | Mouse gut | N/A|
| [Lloyd-Price et al, 2019](https://doi.org/10.1038/s41586-019-1237-9) | Human gut | Ulcerative colitis Vs nonIBD |

