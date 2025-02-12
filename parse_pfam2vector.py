#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 14:59:20 2025

@author: jlu8

script to generate binary pfam vector file from:
    - list of genome pfam results (--genome)
    - lis of sample pfam results (--samples)
    
"""



import re
import os
import glob
from functools import reduce
import argparse
from collections import OrderedDict, Counter
import sys
import warnings
import numpy as np # type: ignore
import pandas as pd # type: ignore
import time
from io import StringIO
import csv



#read pfam clan data file into memory
def read_pfamClanInfo2df(pfam_clan_file):

    pfam_df = pd.read_csv(pfam_clan_file, sep='\t', names = ["Protein_family","class","name","target.name","Description"])

    #rename Protein_family to PfamID
    pfam_df.rename(columns={'Protein_family': 'PfamID'}, inplace=True)

    return pfam_df



## scripts to generate pfam vector from metagenome
def extract_pfams_from_cluster(cluster_file):
    pfam_accessions = []
    
    
    genome_df = pd.read_csv(cluster_file, sep='\t', header=1)
    
    for line in list(genome_df['pfam'].unique()):
        
        if ',' in line:
            line = line.split(",")
            pfam_accessions.extend(line)
        elif line != '-':
            pfam_accessions.append(line)            
        

    return pfam_accessions

## scripts to generate pfam vector from metagenome
def extract_pfams2dict(pfam_clan_df, functional_profiles_folder):
    
    functional_profiles = os.listdir(functional_profiles_folder)
    
    pfam_dict = {pfam: [] for pfam in list(pfam_clan_df['PfamID'].unique())}
    
    obselete_pfams = []
    
    
    genomes_list = [line.split("_")[0] for line in functional_profiles]
    
    genomes_list.insert(0, "PfamID")
    
    for file in functional_profiles:
        
        genome_id = file.split("_")[0]
        cluster_file = os.path.join(functional_profiles_folder, file)
        
        pfam_accessions = extract_pfams_from_cluster(cluster_file)

        #add genome_id to dict        
        for pfam in pfam_accessions:
            if pfam in pfam_dict:
                pfam_dict[pfam].append(genome_id)    
            else:
                obselete_pfams.append(pfam)
                
    obselete_pfams = list(set(obselete_pfams))    

    ## NEED to think about how to deal with deceprecated pfams

    if len(obselete_pfams) > 0:
        print (f"reference metagenome contains a list of {str(len(obselete_pfams))} obselete pfams.")
        

        
    #remove keys with empty lists
    pfam_dict = {k: v for k, v in pfam_dict.items() if v != []}



    return genomes_list, pfam_dict, obselete_pfams
    


## scripts to generate pfam vector from metagenome
def generate_genome_pfam_vector2(pfam_clan_df, genomes_list, pfam_dict, result_genome_matrix_file):
    
    
    # # get list of pfams (pfam-clans) and list of proteins
    pfam_list = list(pfam_clan_df['PfamID'])
    
    with open(result_genome_matrix_file, "w") as rfile:
        title = "\t".join(genomes_list)+"\n"
        rfile.write(title)
        #go through pfam pfam_clan_df and get all genomes
        for i, pfam in enumerate(pfam_list):
            if pfam in pfam_dict:
                rline = ['1' if genome in pfam_dict[pfam] else '0' for genome in genomes_list]
            else:
                rline = ['0' for genome in genomes_list]
                
            rline.insert(0, pfam)
            rline = "\t".join(rline)+"\n"       
            rfile.write(rline)
            

    return i


#generate pfam vectors from folder of pfam interpro files.
def generate_pfam_samples_vector(pfam_clan_df, sample_pfam_vector_folder, pfam_vector_file):
    
    #get list of pfamIDs 
    pfam_clan_ids = list(pfam_clan_df['PfamID'])
    
    #pfam dict - pfam: [samples]
    pfam_dict = {pfam : [] for pfam in pfam_clan_ids}
    
    #assign pfam dict
    pfam_df = pfam_clan_df[['PfamID']]
    
    #get list of all interpro tsv files
    sample_pfam_files = [file  for file in os.listdir(sample_pfam_vector_folder) if file.endswith("pfam.tsv")]
    
    samples = [file.split("_")[0] for file in sample_pfam_files]
    
    title = ["PfamID"]
    title.extend(samples)
    title = "\t".join(title)+"\n"

    obselete_pfams = [] #contains list of obselete pfams
        
    for sample_file in sample_pfam_files:
        assembly_accession = sample_file.split("_")[0]
        sample_file = os.path.join(sample_pfam_vector_folder, sample_file)
        
        #bug? format as csv file
        pfam_df = pd.read_csv(sample_file, sep=',', names=["protein", "PfamID", "pfam_description"])
        
        #get list of pfams from each sample and read into dict
        pfams_list = list(pfam_df["PfamID"])
        
        for pfam in pfams_list:
            # add to dict
            try:
                pfam_dict[pfam].append(assembly_accession)
            except KeyError:
                #how to deal with this?
                obselete_pfams.append(pfam)
                
                # print (f'{pfam} is obselete.') # How do we deal with obselete pfams

    print (f"samples contains a list of {str(len(obselete_pfams))} obselete pfams.")
            
    ## TO DO: How do we deal with obselete proteins
    with open(pfam_vector_file, "w") as pfile:
        pfile.write(title)
        for pfam in pfam_dict:
            if len(pfam_dict[pfam]) != 0:
                rline = ['1' if sample in pfam_dict[pfam] else '0' for sample in samples]
            else:
                rline = ['0' for sample in samples]

            rline.insert(0, pfam)
            rline = "\t".join(rline)+"\n"
            
            pfile.write(rline)
            
        
                
        
    
    
    return obselete_pfams



def main():
    parser = argparse.ArgumentParser(
        description="script creates a binary matrix (presence|absence (1|0) of pfams from either MGnify-analysed assemblies (sample), and/or  MGnify genomes catalogues. The generated vector matrix files can be used directly in MiMiC2 to estimate the microbial genomes that functionally represent a specific biom (synthetic community).",
        usage="python parse_pfam2vector.py --pfam_clan <pfam-clans.tsv> --samples </path/to/folder/of/MGnify-analysed-pfam-anootations> --genomes </path/to/folder/of/MGnify/genomes/catalogue/pan-genome/cluster/results> --prefix <result_prefix>"
        )
    parser.add_argument(
        "--pfam_clan",
        type=str,
        help="full path to Pfam-A.clans.tsv of the pfam database that was used to annotate the selected genomes catalogue and sample(s)",
        required=True
        )

    parser.add_argument(
        "--samples",
        type=str,
        help="fullpath to folder containing pfam results from MGnify-analysed assemblies. Files in folder must have the extension (_FASTA_pfam.tsv)",
        required=False
    )
    
    parser.add_argument(
        "--genomes",
        type=str,
        help="fullpath to folder containing pangeome cluster results from a MGnify genomes catalogue. Files in this folder must be named as MGYGXXXXXXX_clstr.tsv",
        required=False
    )
    
    parser.add_argument(
        "--prefix",
        type=str,
        help="prefix to result matrices",
        required=True,
        default="genome_catalogue_vector"
    )


    args = parser.parse_args()
    
    #assign variables
    pfam_clan_file = args.pfam_clan
    sample_pfam_vector_folder = args.samples
    functional_profiles_folder = args.genomes
    result_prefix = args.prefix    
    
    pfam_clan_df = read_pfamClanInfo2df(pfam_clan_file) #used to create both genome and sample matrices
    
    
    if sample_pfam_vector_folder is not None:
        sample_pfam_vector_file = result_prefix+"_genomes_pfam_vector.tsv"
        samples_obselete_pfams = generate_pfam_samples_vector(pfam_clan_df, sample_pfam_vector_folder, sample_pfam_vector_file)
        print (f"Pfam vector matrix created for {sample_pfam_vector_file}.\n")

    if functional_profiles_folder is not None:
        result_genome_matrix = result_prefix+"_samples_pfam_vector.tsv"
        genomes_list, pfam_dict, genome_obselete_pfams = extract_pfams2dict(pfam_clan_df, functional_profiles_folder)
        generate_genome_pfam_vector2(pfam_clan_df, genomes_list, pfam_dict, result_genome_matrix)
        print (f"Pfam vector matrix created for {result_genome_matrix}.\n")

    if sample_pfam_vector_folder is None and functional_profiles_folder is None:
        #if no folders are specified, return error
        print ("Neither --genomes or --samples folders were specified. Either or both --genomes and --samples must be specified for the functions to work.")

    
if __name__ == "__main__":
    main()