#!/usr/bin/env python3.6
import sys
import glob
import pandas as pd
import operator
import collections
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import itertools
import subprocess
import random
from tqdm import tqdm
import pandas
import matplotlib
from matplotlib import rc
import seaborn as sns
sns.set_palette("hls",15)
matplotlib.rcParams['pdf.fonttype'] = 42
from statistics import mean
import os
import argparse
import warnings
from time import process_time

warnings.filterwarnings("ignore")


location_of_file = os.path.dirname(os.path.realpath(__file__))


location_of_call = os.getcwd()

version = '2025-02-07 (JL)'

# User input
parser = argparse.ArgumentParser(description='MiMiC2 v' + version)

# Options for Sample data
parser.add_argument('-s','--samples', metavar='{INPUT}', required=True, help='Pfam vector file of all metagenomic samples to be studied.')
parser.add_argument('-m','--metadata', metavar='{INPUT}', required=False, help='Metadata file detailing the group assignment of each sample.')
parser.add_argument('--group', metavar='{INPUT}', required=False, help='Name of the group of interest for SynCom creation.')

# Option for Function database
parser.add_argument('-p','--pfam', metavar='{INPUT}',required=True, help='Pfam file e.g. Pfam-A.clans.csv, provided for Pfam v32 in `datasets/core/`')


# Options for Genome database
parser.add_argument('-g','--genomes', metavar='{INPUT}',required=True, help='Pfam vector file of genome collection.')
parser.add_argument('--models', metavar='{INPUT}', required=False, help='Folder containing metabolic models for each genome. Must be provided as RDS files such as those provided by GapSeq. If no folder is provided, than the metabolic modelling step (step 4) will be skipped.')
parser.add_argument('-t','--taxonomy', metavar='{INPUT}',required=False, help='Taxonomic assignment of each genome.')
parser.add_argument('--taxonomiclevel', metavar='{INPUT}',required=False, help='Taxonomic level for filtering (species = s,genus = g, class = c, order = o, phyla = p).')


# Options for SynCom design
parser.add_argument('-c','--consortiasize', metavar='{INT}',required=True,default=10, help='Define the SynCom size the user is after.')
parser.add_argument('--corebias', metavar='{FLOAT}', required=False,default=0.0005, help='The additional weighting provided to Pfams core to the studied group (default = 0.0005).')
parser.add_argument('--groupbias', metavar='{FLOAT}', required=False,default=0.0012, help='The additional weighting provided to Pfams significantly enriched in the studied group (default = 0.0012).')
parser.add_argument('--prevfilt', metavar='{FLOAT}', required=False,default=33.3, help='Prevalence filtering threshold for shortlisting genomes for inclusion in the final SynCom selection (default = 33.3).')


# Options for output
parser.add_argument('-o','--output', metavar='{OUTPUT}', required=True, help='Prefix for all output files e.g. HuSynCom.')

# Advanced user options
parser.add_argument('--iterations', metavar='{INT}',required=False,default=20, help='Change the number of iterations to select sample specific strains in step 1.')
parser.add_argument('--exclusion', metavar='{INPUT}', required=False, help='Provide a list of genome names to be excluded during SynCom selection (in tsv format). Genome names must be the same as the ones included in the genome collection vector file (--genome).')

args, unknown = parser.parse_known_args()

start_time = process_time()

#Variable assignments
consortia_size_wanted = int(args.consortiasize)
reduce_genome_list_number = int(args.iterations)
pfam_file = args.pfam
sample_file = args.samples
genome_file = args.genomes
grouping_file = args.metadata
models_folder = args.models


if args.taxonomy is not None:
	taxonomic_filtering = True
	taxonomic_file = args.taxonomy
	taxonomic_level = args.taxonomiclevel
else:
    taxonomic_filtering = False

	
core_bias = float(args.corebias)
pfam_bias = args.groupbias
prevalence_for_shortlisting = args.prevfilt
group_to_be_studied = args.group # Defines the group we predict for ##### Optional
output_prefix = args.output
			

print ("Reading in pfam vector files. ")
## Read in files
pfams = []
for line in open(pfam_file,'r'):
    timber = line.split('\t')
    pfams.append(timber[0])
print ('Number of pfams studied; ' + str(len(pfams))) #number of pfams in the pfam-clans file used


samples = pd.read_csv(sample_file, sep='\t', index_col='PfamID')

#if exclusion list is provided, genomes (of pathogens) will be exluded from the final vector here.
genomes = pd.read_csv(genome_file, sep='\t', index_col='PfamID')

# Give an option that people can state species to be ignored, but names must be the same as those in the file
if args.exclusion is not None:
	with open(args.exclusion,'r') as file:
		pathogens = [line.rstrip() for line in file]
	genomes = genomes.drop(pathogens, axis=1)
	print ("The following genomes will be excluded from the analysis:", ", ".join(pathogens))
else:
    pathogens = []


# if grouping file is assigned
if grouping_file is not None:
    if grouping_file.endswith("tsv"):
        groups = pd.read_csv(grouping_file, sep='\t', index_col='SampleID')
    else:
        groups = pd.read_csv(grouping_file, sep=',', index_col='SampleID')
    samples_to_be_studied = list(groups.loc[groups['Group'] == group_to_be_studied].index)
else:
    groups = ''
    samples_to_be_studied = list(samples.columns)

print ('Number of genomes studied; ' + str(len(genomes.keys())))

find_models = {}
if models_folder is not None and models_folder != " ":
    for cfile in glob.glob(models_folder + '*.RDS'):
        model_name = cfile.split('/')[-1]
        find_models[model_name] = cfile
    
    print ('Number of GEMs; ' + str(len(find_models.keys())))
else:
      print ('Number of GEMs; ' + str(len(find_models.keys())))
      print (':: Metabolic modelling not assigned. Terminating at Step 4 ::')


print (':: Reading in pfam and vector files')

## Determine sample groups


if isinstance(groups, str):
    print ('Total samples = ', str(len(samples.columns)))           
    print ('All samples will be studied as a single group, without comparative weighting.')
    group_to_be_studied = list(samples.columns) 
    g1_samples = list(samples.columns) 
    g1_bias = []
else:
    g1_bias = []

    g2_bias = []

    g1_samples = []

    g2_samples = []

    group1_is = list(groups['Group'].unique())[0]
    group2_is = list(groups['Group'].unique())[1]
    for SampleID, info in groups.iterrows():
        grouping = info['Group']
        if grouping == group1_is:
            g1_samples.append(SampleID)
        if grouping == group2_is:
            g2_samples.append(SampleID)


    print ('Total samples = ', str(len(samples.columns)))           
    print ('Group 1 = ', group1_is)           
    print ('Number of samples = ', str(len(g1_samples)) + '\n')           

    print ('Group 2 = ', group2_is)           
    print ('Number of samples = ', str(len(g2_samples)))    

    # Set the group to be studied and a consortia predicted for
    group_to_be_used = ''

    if group1_is == group_to_be_studied:
        group_to_be_studied == g1_samples

    if group2_is == group_to_be_studied:
        group_to_be_studied == g2_samples

    for i, j in samples.iterrows():
        g1_pres = 0
        g1_abs = 0
        g2_pres = 0
        g2_abs = 0    

        for sample in g1_samples:
            if j[sample] == 1:
                g1_pres +=1
            else:
                g1_abs +=1

        for sample in g2_samples:
            if j[sample] == 1:
                g2_pres +=1
            else:
                g2_abs +=1

        oddsratio, pvalue = stats.fisher_exact([[g1_pres, g1_abs], [g2_pres, g2_abs]])
        #print ([g1_pres, g1_abs], [g2_pres, g2_abs])
        #print (pvalue)
        if pvalue <0.05:
            g1_ratio = g1_pres / (g1_pres+g1_abs)
            g2_ratio = g2_pres / (g2_pres+g2_abs)
            if g1_ratio > g2_ratio:
                g1_bias.append(j.name)
            else:
                g2_bias.append(j.name)


    print ('Number of Group 1 enriched Pfams; ', len(g1_bias))
    print ('Number of Group 2 enriched Pfams; ',len(g2_bias))

if isinstance(groups, str):
    g1_core_bias = []

    for pfam, samples_info in samples.iterrows():
        g1_pres = 0
        g1_abs = 0


        for sample in g1_samples: # only show info for this group
            if samples_info[sample] == 1:  # If Pfam present in sample
                g1_pres +=1
            else:
                g1_abs +=1

        g1_ratio = g1_pres / (g1_pres+g1_abs)

        if g1_ratio > 0.5:
            g1_core_bias.append(pfam)

    print ('Core proteins in Group 1 = ' , len(g1_core_bias))
    
    
else:    
    g1_core_bias = []

    for pfam, samples_info in samples.iterrows():
        g1_pres = 0
        g1_abs = 0


        for sample in g1_samples: # only show info for this group
            if samples_info[sample] == 1:  # If Pfam present in sample
                g1_pres +=1
            else:
                g1_abs +=1

        g1_ratio = g1_pres / (g1_pres+g1_abs)

        if g1_ratio > 0.5:
            g1_core_bias.append(pfam)

    print ('Core proteins in Group 1 = ' , len(g1_core_bias))

    g2_core_bias = []

    if len (g2_samples) > 0:
        for pfam, samples_info in samples.iterrows():
            g2_pres = 0
            g2_abs = 0


            for sample in g2_samples: # only show info for this group
                if samples_info[sample] == 1:  # If Pfam present in sample
                    g2_pres +=1
                else:
                    g2_abs +=1

            g2_ratio = g2_pres / (g2_pres+g2_abs)

            if g2_ratio > 0.5:
                g2_core_bias.append(pfam)

        print ('Core proteins in Group 2 = ' , len(g2_core_bias))
print ('::: Determined functions weight.')


## Step 1
print (':::: Step 1: Iterative sample-specific selection of strains. ::::')

wanted_pfams = []

iterations_wanted = reduce_genome_list_number

Stored_first_round_consortia = {}
for column in tqdm(samples[samples_to_be_studied].columns): # For each sample        ###########
    sample_data = samples[column]
    sample_specific_genomes = genomes[list((set(genomes.keys())-set(pathogens)))]
    consortia = []
    consortia_matches = []
    consortia_mismatches = []
    cul_accounted = []
    accouted_pfams = []
    accouted_pfams_percentage = 0.0
    samples_pfams = list(sample_data.loc[sample_data > 0].index)
    samples_remaining_pfams = samples_pfams

    #print ('Minimal consortia for sample; ', column)
    #print ('Sample contains this many pfams; ', str(len(samples_pfams)) )

    for iteration in range(0,iterations_wanted): # Create an initial list of 50 consortia members
        genome_scores = {}
        #print ('Genome selection started; ', datetime.now())
        #print ('Number of genomes;', str(len(sample_specific_genomes.columns)))
        for genome in sample_specific_genomes: # Iterate over the remaining genomes each time

            genome_data = sample_specific_genomes[genome] # This section creates the single genome file

            # This is the key difference, here we select to study only those present in the genome
            genome_pfams = list(genome_data.loc[genome_data == 1].index)

            #!!!!!!!!!!!!!!!! This is the time intensive step
            #t1 = [x for x in genome_pfams if x not in accouted_pfams] # Removes those pfams already accounted for
            #genome_pfams = t1

            match_list = list(set(samples_remaining_pfams).intersection(genome_pfams))
            match = 0.0

            # Bias based on core sample functions
            if samples_to_be_studied == g1_samples:
                bias_core_score = len(list(set(g1_core_bias).intersection(genome_pfams))) * core_bias
                bias_spef_score = len(list(set(g1_bias).intersection(genome_pfams))) * pfam_bias
            elif samples_to_be_studied == g2_samples:
                bias_core_score = len(list(set(g2_core_bias).intersection(genome_pfams))) * core_bias
                bias_spef_score = len(list(set(g2_bias).intersection(genome_pfams))) * pfam_bias
            # Bias based on core genome functions
            for i in match_list:
                #match += pfams_prev[i]
                match +=1
            #print (match)
            mismatch = len(list(genome_data.loc[genome_data == 1].index))-len(match_list)
            #print (mismatch)

            try:
                score = (match/ (match+mismatch)) + bias_core_score + bias_spef_score
                genome_scores[genome] = [score, match, mismatch]
            except:
                print ("FAILED to be studied.")
                print (column, len(accouted_pfams), len(samples_remaining_pfams))
                print (genome, len(genome_pfams),match, mismatch)

        #print ('Genome selection ended; ', datetime.now())
        max_match = max(genome_scores, key=genome_scores.get) # Identify the best consortia member this round
        consortia.append(max_match) # Add to list
        #print ('Added consortia member; ', max_match)

        # At this point we have selected the consortia member so the remainder readies the datasets for the next round

        # Obtain best matches data
        genome_data = sample_specific_genomes[max_match] # This section creates the single genome file
        genome_pfams = list(genome_data.loc[genome_data == 1].index)
        score, match, mismatch = genome_scores[max_match]
        #print (match, mismatch, score)

        # Calculate cul. accounted pfams %
        for i in list(set(samples_remaining_pfams).intersection(genome_pfams)):
            accouted_pfams.append(i)
        try:
            accouted_pfams_percentage = (len(accouted_pfams)/ float(len(samples_pfams)))*100
        except ZeroDivisionError:
            if len(accouted_pfams) == 0:
                if len(samples_pfams) == 0:
                    accouted_pfams_percentage = 0.0
            
        #print ('Culumalative pfams accounted for; ', str(accouted_pfams_percentage))

        # Remove included genome from genome list
        t1 = sample_specific_genomes.drop(list(set(samples_remaining_pfams).intersection(genome_pfams)))
        t2 = t1.drop(columns=[max_match])
        sample_specific_genomes = t2

        # This section removes pfams already counted for and the selected taxa
        s1 = [x for x in samples_remaining_pfams if x not in list(set(samples_pfams).intersection(genome_pfams))]
        samples_remaining_pfams = s1
        #print ("Unaccounted for pfams; ", str(len(samples_remaining_pfams)))

        # Store statistics
        consortia_matches.append(match)
        consortia_mismatches.append(mismatch)
        cul_accounted.append(accouted_pfams_percentage)
        #print ('Selection made for round; ', str(iteration))

    cul_mismatches = []
    tot_mismatch = 0
    for i in consortia_mismatches:
        tot_mismatch += i
        cul_mismatches.append(tot_mismatch)

    #print(round(kneedle.knee, 3))
    #print(round(kneedle.elbow, 3))
    #kneedle.plot_knee_normalized()
    #kneedle.plot_knee()
    Stored_first_round_consortia[column] = [consortia_matches, consortia_mismatches, cul_accounted,cul_mismatches, consortia]

    #fig = plt.figure()
    #ax = plt.axes()
    #ax.plot(consortia_mismatches)


    #fig = plt.figure()
    #ax = plt.axes()
    #ax.plot(cul_mismatches)

print (':::: Step 1 complete: Iterative sample-specific selection of strains.')

## Step 2
print ('::::: Step 2: Strain reduction and consortia filtering.')


genome_prev = {}

for sample, data in Stored_first_round_consortia.items():
    sample_Db = []
    #print (line.split('\t'))
    min_size = data[4]
    numbs = 0
    for gen in data[4]:
        numbs +=1
        #print (numbs)
        if numbs <= reduce_genome_list_number: 
            #print (numbs)
            if gen in genome_prev:
                genome_prev[gen] +=1
            else:
                genome_prev[gen] = 1

                    
                    
                    
print ('Total space of seleted species;' + str(len(genome_prev)))

#### Determine number of samples required to be positive in
if samples_to_be_studied == g1_samples:
    samples_needed_in = float(float(prevalence_for_shortlisting)/100)*len(g1_samples)
elif samples_to_be_studied == g2_samples:
    samples_needed_in = float(float(prevalence_for_shortlisting)/100)*len(g2_samples)

wanted_prev = []


for k,v in genome_prev.items():
    if v >= samples_needed_in:                              # Filter samples based on percentage stated
        wanted_prev.append(k)

        
print ('Prevalence filtered species;' + str(len(wanted_prev)))
wanted = list((set(wanted_prev))-set(pathogens))
wanted_prev = wanted


outputting_filt = location_of_call + '/'+ output_prefix + '-Studied_strains.tsv'

with open(outputting_filt, 'w') as file:
    file.write("Species\n")
#    for species in wanted_prev:
#        file.write(species + '\n')

consortia_size = consortia_size_wanted


if len(wanted_prev) < consortia_size:
    print ('ERROR: Insufficient strains have passed filtering, change the prevalence filter to retain sufficient diversity for consortia creation.')
    sys.exit()

all_consortia_combinations = itertools.combinations(genomes[wanted_prev].keys(), consortia_size)    ###############


consortia_scores = {}
for consortia in all_consortia_combinations:
    consortia_scores[consortia] = []

print ('Number of combinations; ', str(len(consortia_scores.keys())))
if taxonomic_filtering == True:
    print ('Filtered based on taxonomic level; ', taxonomic_level)
    filtered_consortia_scores = {}
    taxonomy = pd.read_csv(taxonomic_file,sep='\t',index_col='user_genome')
    for comb in consortia_scores:
        taxonomies = []
        #print (comb)
        for i in comb:
            d, p, c, o, f, g, s = taxonomy.loc[i]['classification'].split(';')
            #print (i, s)
            if taxonomic_level == 's':
                taxonomies.append(s)
            if taxonomic_level == 'g':
                taxonomies.append(g)
            if taxonomic_level == 'f':
                taxonomies.append(f)
            if taxonomic_level == 'o':
                taxonomies.append(o)
            if taxonomic_level == 'c':
                taxonomies.append(c)
            if taxonomic_level == 'p':
                taxonomies.append(p)
        #print (taxonomies)
        #print (len(set(taxonomies)))
        if len(set(taxonomies)) == len(comb):          # 
            filtered_consortia_scores[tuple(comb)] = []

consortia_scores = filtered_consortia_scores
print ('Number of combinations after filtering; ', len(consortia_scores.keys()))

print ('::::: Step 2 complete: Strain reduction and consortia filtering.')

## Step 3
print (':::::: Step 3: Group scoring of consortia.')
for column in tqdm(list(samples[samples_to_be_studied].columns)): # For each sample        ###########
    consortia_studied = 0
    redundancies_removed = 0
    sample_data = samples[column]
    matches = []
    mismatches = []
    consortia = []
    consortia_matches = []
    consortia_mismatches = []
    cul_accounted = []
    accouted_pfams = []
    accouted_pfams_percentage = 0.0
    samples_pfams = list(sample_data.loc[sample_data == 1].index)
    samples_remaining_pfams = samples_pfams
    #print ('Total number of Pfams; ', len(samples_pfams), ' Remaining Pfams; ', len(samples_remaining_pfams))    
    sample_specific_genomes = genomes

    best_consortia = {}
    failed_consortia = []

    for combination in consortia_scores.keys():       ###########
        # Prevent redundant analysis!!!!!!!!
        combinations_pfams = []
        for genome in combination:   
            genome_data = genomes[genome] # This section creates the single genome file
            redundant_functional_list = combinations_pfams + list(genome_data.loc[genome_data == 1].index)
            combinations_pfams = list(set(redundant_functional_list))


        match = 0
        mismatch = 0

        match = len(list(set(samples_remaining_pfams).intersection(combinations_pfams))) # Identifies the common 
        #for i in combinations_pfams:
            #match += pfams_prev[i]
            
        # Bias based on core sample functions
        if samples_to_be_studied == g1_samples:
            bias_core_score = len(list(set(g1_core_bias).intersection(combinations_pfams))) * core_bias
            bias_spef_score = len(list(set(g1_bias).intersection(combinations_pfams))) * pfam_bias
        elif samples_to_be_studied == g2_samples:
            bias_core_score = len(list(set(g2_core_bias).intersection(combinations_pfams))) * core_bias
            bias_spef_score = len(list(set(g2_bias).intersection(combinations_pfams))) * pfam_bias
                
        mismatch = len(combinations_pfams)-match

        score = (match/ (match+mismatch)) + bias_core_score + bias_spef_score
        prev = consortia_scores[combination]
        prev.append(score)
        
        consortia_scores[combination] = prev

consortia_mean_Scores = {}
for k,v in consortia_scores.items():
    consortia_mean_Scores[k] = np.mean(v)
max_match = max(consortia_mean_Scores, key=consortia_mean_Scores.get) # Identify the best consortia member this round

print (":::::: Consortia Summary ::::::")
print ('Best consortia is; ' + str(max_match))
print ('Consortia-wide score; ', consortia_mean_Scores[max_match])
print ('Consortia tested; ' + str(len(consortia_mean_Scores.keys())))


consortia_selected = max_match               ##############

outputting_all_Scores = location_of_call + '/'+ output_prefix + '-Consortia_scores.tsv'

with open(outputting_all_Scores, 'w') as file:
    file.write('SynCom_members\tGroup_score\n')
    for k,v in consortia_mean_Scores.items():
        file.write(str(k) + '\t' + str(v) + '\n')

end_time = process_time()

step4_start_time = process_time()

print (':::::: Step 3 complete: Group scoring of consortia.')



# Step 4
print ('::::::: Step 4: Metabolic modelling and interaction analysis.')
# Need to write a for loop so only when these two flags are True, does it end the code
# This week mean iterations occur until a stable consortia is found
passed_phase_one_stability = False
passed_phase_two_stability = False

#output files
community_outputfile = location_of_call + '/'+ output_prefix + '-Community-wide_response.tsv'
paired_outputfile = location_of_call + '/'+ output_prefix +'-Paired_Interactions.tsv'

paired_title = "\t".join(["taxa1 and taxa2", "interaction"])+"\n"
community_title = "\t".join(["species", "consortia_vs_single_growth"])+"\n"

modelling_iteration = 0

individual_growths = {}

if models_folder is not None and models_folder != " ":
    while passed_phase_two_stability == False:
        #max_match = max(consortia_mean_Scores, key=consortia_mean_Scores.get) # Identify the best consortia member this round
        #consortia_selected = max_match               ##############

        models = [] # list of models for the selected consortia

        for species in consortia_selected:
            models.append(find_models[species + '.RDS'])

        all_data_stored = {} # Store consortia predictions


        if consortia_size_wanted == 2:
            for listi,individual in enumerate(models):
                tmp_models = list(models)
                tmp_models.pop(listi)

                name1 = models[listi].split('/')[-1:][0].split('.RDS')[0]
                name2 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]

                print ('Currently growing; ', name1)


                # !!!!!!!!!!!!! Here, the species being studied is given first, then the other ones, then all the names
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size2/Single_Growth.R', location_of_call + '/' , models[listi], tmp_models[0], name1,name2])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                growth_curve = []
                ln_num = 0
                for line in open(location_of_call + '/' + 'statistics.csv', 'r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1]
                        growth_curve.append(int(line.replace('\n','').split(',')[3]))
                individual_growths[name1] = growth_curve

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv '])

            print ('Individual growth curves generated for ;' + str(len(individual_growths)))

                
            all_consortia_combinations = list(itertools.combinations(models, 2))

            Paired_growths = {}

            for comb in all_consortia_combinations:
                tmp_models = list(models[:])
                tmp_models.remove(list(comb)[0])
                tmp_models.remove(list(comb)[1])

                name1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                name2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                print ('Strains being grown together;')
                print( list(comb)[0].replace('"','').split('/')[-1:][0].split('.RDS')[0], list(comb)[1].replace('"','').split('/')[-1:][0].split('.RDS')[0], '; names=(', name1,',',name2,')') # Only highlighting those species being paired
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size2/Paired_Growth.R', location_of_call + '/' , list(comb)[0], list(comb)[1], name1,name2])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                ln_num = 0
                growth_info = {}
                growth_info[list(comb)[0].split('/')[-1:][0].split('.RDS')[0]] = []
                growth_info[list(comb)[1].split('/')[-1:][0].split('.RDS')[0]] = []
                for line in open(location_of_call + '/' + 'statistics.csv', 'r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1].replace('"','').split('/')[-1:][0].split('.RDS')[0]
                        growth = int(line.replace('\n','').split(',')[3])
                        growth_info[species].append(growth) 
                Paired_growths[comb] = growth_info

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv '])

            ### Identify paired interactions
            with open(paired_outputfile, 'w') as outputting:
                outputting.write(paired_title)
                for comb, data in Paired_growths.items():
                    print ('Calculating interaction between; ')
                    #print (comb)
                    species_1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                    species_2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                    print (species_1)
                    print (species_2)

                    # Calculate interaction type for species1

                    single1 = individual_growths[species_1][-1:][0]
                    paired1 = data[species_1][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    if perc1 >10.0:
                        interaction1 = 'positive'
                    elif perc1 < -10.0:
                        interaction1 = 'negative'
                    else:
                        interaction1 = 'none'


                    # Calculate interaction type for species2

                    single2 = individual_growths[species_2][-1:][0]
                    paired2 = data[species_2][-1:][0]

                    perc2 = ((single2/paired2)*100)-100
                    print (perc1, perc2)

                    if perc2 >10.0:
                        interaction2 = 'positive'
                    elif perc2 < -10.0:
                        interaction2 = 'negative'
                    else:
                        interaction2 = 'none'


                    # Overall interaction type

                    if interaction1 == 'positive':
                        if interaction2 == 'positive':
                            overall_interaction = 'mutualism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'none':
                            overall_interaction = 'commensalism'

                    elif interaction1 == 'negative':
                        if interaction2 == 'positive':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'competition'
                        elif interaction2 == 'none':
                            overall_interaction = 'amensalism'


                    elif interaction1 == 'none':
                        if interaction2 == 'positive':
                            overall_interaction = 'commensalism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'amensalism'
                        elif interaction2 == 'none':
                            overall_interaction = 'neutralism'

                    outputting.write(species_1 + ' and ' + species_2 + '\t ' + overall_interaction + '\n')
                    print (species_1 + ' and ' + species_2 + '; ' + overall_interaction + '\n')

            print ('Running whole community model.')
            Community_growth = {}

            name1 = models[0].split('/')[-1:][0].split('.RDS')[0]
            name2 = models[1].split('/')[-1:][0].split('.RDS')[0]
            try:
                subprocess.run(['rm', location_of_call + '/' + 'statistics.csv'])
            except FileNotFoundError:
                lol = 0

            #print( list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1], name1,name2,name3,name4)
            subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size2/Combined_Growth.R', location_of_call + '/' , models[0], models[1], name1,name2])
            # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4

            ln_num = 0
            growth_info = {}

            Community_growth[models[0].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[1].split('/')[-1:][0].split('.RDS')[0]] = []

            for line in open(location_of_call + '/' + 'statistics.csv', 'r'): #change
                ln_num +=1
                if ln_num >1:
                    print (",".join(line.split(',')[:-1][1:]))
                    species = line.split(',')[1].replace('"','')
                    growth = int(line.replace('\n','').split(',')[3])
                    Community_growth[species].append(growth) 

            #subprocess.run(['rm ' + location_of_call + '/' + 'statistics.csv'])

            with open(community_outputfile, 'w') as outputting:
                outputting.write(community_title)

                species_1 = models[0].split('/')[-1:][0].split('.RDS')[0]
                species_2 = models[1].split('/')[-1:][0].split('.RDS')[0]


                for species in [species_1,species_2]:
                    single1 = individual_growths[species][-1:][0]
                    paired1 = Community_growth[species][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    outputting.write(species + '\t' + str(perc1) + '%\n')


            ## Store the data

            all_data_stored[tuple(models)] = [individual_growths,Paired_growths,Community_growth]

            #print ('Completed consortia; ', str(len(all_data_stored.keys())))


        # ## Three species

        # In[45]:


        if consortia_size_wanted == 3:
            for listi,individual in enumerate(models):
                tmp_models = list(models)
                tmp_models.pop(listi)

                name1 = models[listi].split('/')[-1:][0].split('.RDS')[0]
                name2 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]

                print ('Currently growing; ', name1)


                # !!!!!!!!!!!!! Here, the species being studied is given first, then the other ones, then all the names
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size3/Single_Growth.R', location_of_call + '/' , models[listi], tmp_models[0],tmp_models[1], name1,name2,name3])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                growth_curve = []
                ln_num = 0
                for line in open(location_of_call + '/' + 'statistics.csv', 'r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1]
                        growth_curve.append(int(line.replace('\n','').split(',')[3]))
                individual_growths[name1] = growth_curve

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv '])

            print ('Individual growth curves generated for ;' + str(len(individual_growths)))

                
            all_consortia_combinations = list(itertools.combinations(models, 2))

            Paired_growths = {}

            for comb in all_consortia_combinations:
                tmp_models = list(models[:])
                tmp_models.remove(list(comb)[0])
                tmp_models.remove(list(comb)[1])

                name1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                name2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                print ('Strains being grown together;')
                print( list(comb)[0].replace('"','').split('/')[-1:][0].split('.RDS')[0], list(comb)[1].replace('"','').split('/')[-1:][0].split('.RDS')[0], '; names=(', name1,',',name2,')') # Only highlighting those species being paired
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size3/Paired_Growth.R', location_of_call + '/' , list(comb)[0], list(comb)[1],tmp_models[0], name1,name2,name3])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                ln_num = 0
                growth_info = {}
                growth_info[list(comb)[0].split('/')[-1:][0].split('.RDS')[0]] = []
                growth_info[list(comb)[1].split('/')[-1:][0].split('.RDS')[0]] = []
                for line in open(location_of_call + '/' + 'statistics.csv', 'r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1].replace('"','').split('/')[-1:][0].split('.RDS')[0]
                        growth = int(line.replace('\n','').split(',')[3])
                        growth_info[species].append(growth) 
                Paired_growths[comb] = growth_info

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv '])

            

            ### Identify paired interactions
            with open(paired_outputfile, 'w') as outputting:
                outputting.write(paired_title)
                for comb, data in Paired_growths.items():
                    print ('Calculating interaction between; ')
                    #print (comb)
                    species_1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                    species_2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                    print (species_1)
                    print (species_2)

                    # Calculate interaction type for species1

                    single1 = individual_growths[species_1][-1:][0]
                    paired1 = data[species_1][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    if perc1 >10.0:
                        interaction1 = 'positive'
                    elif perc1 < -10.0:
                        interaction1 = 'negative'
                    else:
                        interaction1 = 'none'


                    # Calculate interaction type for species2

                    single2 = individual_growths[species_2][-1:][0]
                    paired2 = data[species_2][-1:][0]

                    perc2 = ((single2/paired2)*100)-100
                    print (perc1, perc2)

                    if perc2 >10.0:
                        interaction2 = 'positive'
                    elif perc2 < -10.0:
                        interaction2 = 'negative'
                    else:
                        interaction2 = 'none'


                    # Overall interaction type

                    if interaction1 == 'positive':
                        if interaction2 == 'positive':
                            overall_interaction = 'mutualism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'none':
                            overall_interaction = 'commensalism'

                    elif interaction1 == 'negative':
                        if interaction2 == 'positive':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'competition'
                        elif interaction2 == 'none':
                            overall_interaction = 'amensalism'


                    elif interaction1 == 'none':
                        if interaction2 == 'positive':
                            overall_interaction = 'commensalism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'amensalism'
                        elif interaction2 == 'none':
                            overall_interaction = 'neutralism'

                    outputting.write(species_1 + ' and ' + species_2 + '\t ' + overall_interaction + '\n')
                    print (species_1 + ' and ' + species_2 + '; ' + overall_interaction + '\n')

            
            print ('Running whole community model.')
            Community_growth = {}

            name1 = models[0].split('/')[-1:][0].split('.RDS')[0]
            name2 = models[1].split('/')[-1:][0].split('.RDS')[0]
            name3 = models[2].split('/')[-1:][0].split('.RDS')[0]
            try:
                subprocess.run(['rm', location_of_call + '/' + 'statistics.csv'])
            except FileNotFoundError:
                lol = 0

            #print( list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1], name1,name2,name3,name4)
            subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size3/Combined_Growth.R', location_of_call + '/' , models[0], models[1],models[2], name1,name2,name3])
            # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4

            ln_num = 0
            growth_info = {}

            Community_growth[models[0].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[1].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[2].split('/')[-1:][0].split('.RDS')[0]] = []

            for line in open(location_of_call + '/' + 'statistics.csv''r'): #change
                ln_num +=1
                if ln_num >1:
                    print (",".join(line.split(',')[:-1][1:]))
                    species = line.split(',')[1].replace('"','')
                    growth = int(line.replace('\n','').split(',')[3])
                    Community_growth[species].append(growth) 

            #subprocess.run(['rm ' + location_of_call + '/' + 'statistics.csv'])

            with open(community_outputfile, 'w') as outputting:
                outputting.write(community_title)

                species_1 = models[0].split('/')[-1:][0].split('.RDS')[0]
                species_2 = models[1].split('/')[-1:][0].split('.RDS')[0]
                species_3 = models[2].split('/')[-1:][0].split('.RDS')[0]


                for species in [species_1,species_2,species_3]:
                    single1 = individual_growths[species][-1:][0]
                    paired1 = Community_growth[species][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    outputting.write(species + '\t' + str(perc1) + '%\n')


            ## Store the data

            all_data_stored[tuple(models)] = [individual_growths,Paired_growths,Community_growth]

            #print ('Completed consortia; ', str(len(all_data_stored.keys())))


        # ## Four species

        # In[ ]:


        if consortia_size_wanted == 4:
            for listi,individual in enumerate(models):
                tmp_models = list(models)
                tmp_models.pop(listi)

                name1 = models[listi].split('/')[-1:][0].split('.RDS')[0]
                name2 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]

                print ('Currently growing; ', name1)


                # !!!!!!!!!!!!! Here, the species being studied is given first, then the other ones, then all the names
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size4/Single_Growth.R', location_of_call + '/' , models[listi], tmp_models[0],tmp_models[1],tmp_models[2], name1,name2,name3,name4])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                growth_curve = []
                ln_num = 0
                for line in open(location_of_call + '/' + 'statistics.csv', 'r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1]
                        growth_curve.append(int(line.replace('\n','').split(',')[3]))
                individual_growths[name1] = growth_curve

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv '])

            print ('Individual growth curves generated for ;' + str(len(individual_growths)))

                
            all_consortia_combinations = list(itertools.combinations(models, 2))

            Paired_growths = {}

            for comb in all_consortia_combinations:
                tmp_models = list(models[:])
                tmp_models.remove(list(comb)[0])
                tmp_models.remove(list(comb)[1])

                name1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                name2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                print ('Strains being grown together;')
                print( list(comb)[0].replace('"','').split('/')[-1:][0].split('.RDS')[0], list(comb)[1].replace('"','').split('/')[-1:][0].split('.RDS')[0], '; names=(', name1,',',name2,')') # Only highlighting those species being paired
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size4/Paired_Growth.R', location_of_call + '/' , list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1], name1,name2,name3,name4])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                ln_num = 0
                growth_info = {}
                growth_info[list(comb)[0].split('/')[-1:][0].split('.RDS')[0]] = []
                growth_info[list(comb)[1].split('/')[-1:][0].split('.RDS')[0]] = []
                for line in open(location_of_call + '/' + 'statistics.csv', 'r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1].replace('"','').split('/')[-1:][0].split('.RDS')[0]
                        growth = int(line.replace('\n','').split(',')[3])
                        growth_info[species].append(growth) 
                Paired_growths[comb] = growth_info

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv '])

            ### Identify paired interactions
            with open(paired_outputfile, 'w') as outputting:
                outputting.write(paired_title)
                for comb, data in Paired_growths.items():
                    print ('Calculating interaction between; ')
                    #print (comb)
                    species_1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                    species_2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                    print (species_1)
                    print (species_2)

                    # Calculate interaction type for species1

                    single1 = individual_growths[species_1][-1:][0]
                    paired1 = data[species_1][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    if perc1 >10.0:
                        interaction1 = 'positive'
                    elif perc1 < -10.0:
                        interaction1 = 'negative'
                    else:
                        interaction1 = 'none'


                    # Calculate interaction type for species2

                    single2 = individual_growths[species_2][-1:][0]
                    paired2 = data[species_2][-1:][0]

                    perc2 = ((single2/paired2)*100)-100
                    print (perc1, perc2)

                    if perc2 >10.0:
                        interaction2 = 'positive'
                    elif perc2 < -10.0:
                        interaction2 = 'negative'
                    else:
                        interaction2 = 'none'


                    # Overall interaction type

                    if interaction1 == 'positive':
                        if interaction2 == 'positive':
                            overall_interaction = 'mutualism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'none':
                            overall_interaction = 'commensalism'

                    elif interaction1 == 'negative':
                        if interaction2 == 'positive':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'competition'
                        elif interaction2 == 'none':
                            overall_interaction = 'amensalism'


                    elif interaction1 == 'none':
                        if interaction2 == 'positive':
                            overall_interaction = 'commensalism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'amensalism'
                        elif interaction2 == 'none':
                            overall_interaction = 'neutralism'

                    outputting.write(species_1 + ' and ' + species_2 + '\t ' + overall_interaction + '\n')
                    print (species_1 + ' and ' + species_2 + '; ' + overall_interaction + '\n')

            
            print ('Running whole community model.')
            Community_growth = {}

            name1 = models[0].split('/')[-1:][0].split('.RDS')[0]
            name2 = models[1].split('/')[-1:][0].split('.RDS')[0]
            name3 = models[2].split('/')[-1:][0].split('.RDS')[0]
            name4 = models[3].split('/')[-1:][0].split('.RDS')[0]
            try:
                subprocess.run(['rm', location_of_call + '/' + 'statistics.csv'])
            except FileNotFoundError:
                lol = 0

            #print( list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1], name1,name2,name3,name4)
            subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size4/Combined_Growth.R', location_of_call + '/' , models[0], models[1],models[2],models[3], name1,name2,name3,name4])
            # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4

            ln_num = 0
            growth_info = {}

            Community_growth[models[0].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[1].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[2].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[3].split('/')[-1:][0].split('.RDS')[0]] = []

            for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                ln_num +=1
                if ln_num >1:
                    print (",".join(line.split(',')[:-1][1:]))
                    species = line.split(',')[1].replace('"','')
                    growth = int(line.replace('\n','').split(',')[3])
                    Community_growth[species].append(growth) 

            #subprocess.run(['rm ' + location_of_call + '/' + 'statistics.csv'])
            with open(community_outputfile, 'w') as outputting:
                outputting.write(community_title)

                species_1 = models[0].split('/')[-1:][0].split('.RDS')[0]
                species_2 = models[1].split('/')[-1:][0].split('.RDS')[0]
                species_3 = models[2].split('/')[-1:][0].split('.RDS')[0]
                species_4 = models[3].split('/')[-1:][0].split('.RDS')[0]


                for species in [species_1,species_2,species_3,species_4]:
                    single1 = individual_growths[species][-1:][0]
                    paired1 = Community_growth[species][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    outputting.write(species + '\t' + str(perc1) + '%\n')


            ## Store the data

            all_data_stored[tuple(models)] = [individual_growths,Paired_growths,Community_growth]

            #print ('Completed consortia; ', str(len(all_data_stored.keys())))


        # # Five species

        # In[112]:


        if consortia_size_wanted == 5:
            for listi,individual in enumerate(models):
                tmp_models = list(models)
                tmp_models.pop(listi)

                name1 = models[listi].split('/')[-1:][0].split('.RDS')[0]
                name2 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]

                print ('Currently growing; ', name1)


                # !!!!!!!!!!!!! Here, the species being studied is given first, then the other ones, then all the names
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size5/Single_Growth.R', location_of_call + '/' , models[listi], tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3], name1,name2,name3,name4,name5])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                growth_curve = []
                ln_num = 0
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1]
                        growth_curve.append(int(line.replace('\n','').split(',')[3]))
                individual_growths[name1] = growth_curve

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv '])

            print ('Individual growth curves generated for ;' + str(len(individual_growths)))

                
            all_consortia_combinations = list(itertools.combinations(models, 2))

            Paired_growths = {}

            for comb in all_consortia_combinations:
                tmp_models = list(models[:])
                tmp_models.remove(list(comb)[0])
                tmp_models.remove(list(comb)[1])

                name1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                name2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                print ('Strains being grown together;')
                print( list(comb)[0].replace('"','').split('/')[-1:][0].split('.RDS')[0], list(comb)[1].replace('"','').split('/')[-1:][0].split('.RDS')[0], '; names=(', name1,',',name2,')') # Only highlighting those species being paired
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size5/Paired_Growth.R', location_of_call + '/' , list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1],tmp_models[2], name1,name2,name3,name4,name5])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                ln_num = 0
                growth_info = {}
                growth_info[list(comb)[0].split('/')[-1:][0].split('.RDS')[0]] = []
                growth_info[list(comb)[1].split('/')[-1:][0].split('.RDS')[0]] = []
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1].replace('"','').split('/')[-1:][0].split('.RDS')[0]
                        growth = int(line.replace('\n','').split(',')[3])
                        growth_info[species].append(growth) 
                Paired_growths[comb] = growth_info

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv '])

            ### Identify paired interactions
            with open(paired_outputfile, 'w') as outputting:
                outputting.write(paired_title)
                for comb, data in Paired_growths.items():
                    print ('Calculating interaction between; ')
                    #print (comb)
                    species_1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                    species_2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                    print (species_1)
                    print (species_2)

                    # Calculate interaction type for species1

                    single1 = individual_growths[species_1][-1:][0]
                    paired1 = data[species_1][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    if perc1 >10.0:
                        interaction1 = 'positive'
                    elif perc1 < -10.0:
                        interaction1 = 'negative'
                    else:
                        interaction1 = 'none'


                    # Calculate interaction type for species2

                    single2 = individual_growths[species_2][-1:][0]
                    paired2 = data[species_2][-1:][0]

                    perc2 = ((single2/paired2)*100)-100
                    print (perc1, perc2)

                    if perc2 >10.0:
                        interaction2 = 'positive'
                    elif perc2 < -10.0:
                        interaction2 = 'negative'
                    else:
                        interaction2 = 'none'


                    # Overall interaction type

                    if interaction1 == 'positive':
                        if interaction2 == 'positive':
                            overall_interaction = 'mutualism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'none':
                            overall_interaction = 'commensalism'

                    elif interaction1 == 'negative':
                        if interaction2 == 'positive':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'competition'
                        elif interaction2 == 'none':
                            overall_interaction = 'amensalism'


                    elif interaction1 == 'none':
                        if interaction2 == 'positive':
                            overall_interaction = 'commensalism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'amensalism'
                        elif interaction2 == 'none':
                            overall_interaction = 'neutralism'

                    outputting.write(species_1 + ' and ' + species_2 + '\t ' + overall_interaction + '\n')
                    print (species_1 + ' and ' + species_2 + '; ' + overall_interaction + '\n')

            outputting.close()

            
            print ('Running whole community model.')
            Community_growth = {}

            name1 = models[0].split('/')[-1:][0].split('.RDS')[0]
            name2 = models[1].split('/')[-1:][0].split('.RDS')[0]
            name3 = models[2].split('/')[-1:][0].split('.RDS')[0]
            name4 = models[3].split('/')[-1:][0].split('.RDS')[0]
            name5 = models[4].split('/')[-1:][0].split('.RDS')[0]
            try:
                subprocess.run(['rm', location_of_call + '/' + 'statistics.csv'])
            except FileNotFoundError:
                lol = 0

            #print( list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1], name1,name2,name3,name4)
            subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size5/Combined_Growth.R', location_of_call + '/' , models[0], models[1],models[2],models[3],models[4], name1,name2,name3,name4,name5])
            # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4

            ln_num = 0
            growth_info = {}

            Community_growth[models[0].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[1].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[2].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[3].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[4].split('/')[-1:][0].split('.RDS')[0]] = []

            for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                ln_num +=1
                if ln_num >1:
                    print (",".join(line.split(',')[:-1][1:]))
                    species = line.split(',')[1].replace('"','')
                    growth = int(line.replace('\n','').split(',')[3])
                    Community_growth[species].append(growth) 

            #subprocess.run(['rm ' + location_of_call + '/' + 'statistics.csv'])
            with open(community_outputfile, 'w') as outputting:
                outputting.write(community_title)

                species_1 = models[0].split('/')[-1:][0].split('.RDS')[0]
                species_2 = models[1].split('/')[-1:][0].split('.RDS')[0]
                species_3 = models[2].split('/')[-1:][0].split('.RDS')[0]
                species_4 = models[3].split('/')[-1:][0].split('.RDS')[0]
                species_5 = models[4].split('/')[-1:][0].split('.RDS')[0]


                for species in [species_1,species_2,species_3,species_4,species_5]:
                    single1 = individual_growths[species][-1:][0]
                    paired1 = Community_growth[species][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    outputting.write(species + '\t' + str(perc1) + '%\n')



            ## Store the data

            all_data_stored[tuple(models)] = [individual_growths,Paired_growths,Community_growth]

            #print ('Completed consortia; ', str(len(all_data_stored.keys())))


        # # Six species

        # In[ ]:


        if consortia_size_wanted == 6:
            for listi,individual in enumerate(models):
                tmp_models = list(models)
                tmp_models.pop(listi)

                name1 = models[listi].split('/')[-1:][0].split('.RDS')[0]
                name2 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]

                print ('Currently growing; ', name1)


                # !!!!!!!!!!!!! Here, the species being studied is given first, then the other ones, then all the names
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size6/Single_Growth.R', location_of_call + '/' , models[listi], tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4], name1,name2,name3,name4,name5,name6])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                growth_curve = []
                ln_num = 0
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1]
                        growth_curve.append(int(line.replace('\n','').split(',')[3]))
                individual_growths[name1] = growth_curve

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv '])

            print ('Individual growth curves generated for ;' + str(len(individual_growths)))

                
            all_consortia_combinations = list(itertools.combinations(models, 2))

            Paired_growths = {}

            for comb in all_consortia_combinations:
                tmp_models = list(models[:])
                tmp_models.remove(list(comb)[0])
                tmp_models.remove(list(comb)[1])

                name1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                name2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                print ('Strains being grown together;')
                print( list(comb)[0].replace('"','').split('/')[-1:][0].split('.RDS')[0], list(comb)[1].replace('"','').split('/')[-1:][0].split('.RDS')[0], '; names=(', name1,',',name2,')') # Only highlighting those species being paired
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size6/Paired_Growth.R', location_of_call + '/' , list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3], name1,name2,name3,name4,name5,name6])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                ln_num = 0
                growth_info = {}
                growth_info[list(comb)[0].split('/')[-1:][0].split('.RDS')[0]] = []
                growth_info[list(comb)[1].split('/')[-1:][0].split('.RDS')[0]] = []
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1].replace('"','').split('/')[-1:][0].split('.RDS')[0]
                        growth = int(line.replace('\n','').split(',')[3])
                        growth_info[species].append(growth) 
                Paired_growths[comb] = growth_info

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv '])

            ### Identify paired interactions
            with open(paired_outputfile, 'w') as outputting:
                outputting.write(paired_title)
                for comb, data in Paired_growths.items():
                    print ('Calculating interaction between; ')
                    #print (comb)
                    species_1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                    species_2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                    print (species_1)
                    print (species_2)

                    # Calculate interaction type for species1

                    single1 = individual_growths[species_1][-1:][0]
                    paired1 = data[species_1][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    if perc1 >10.0:
                        interaction1 = 'positive'
                    elif perc1 < -10.0:
                        interaction1 = 'negative'
                    else:
                        interaction1 = 'none'


                    # Calculate interaction type for species2

                    single2 = individual_growths[species_2][-1:][0]
                    paired2 = data[species_2][-1:][0]

                    perc2 = ((single2/paired2)*100)-100
                    print (perc1, perc2)

                    if perc2 >10.0:
                        interaction2 = 'positive'
                    elif perc2 < -10.0:
                        interaction2 = 'negative'
                    else:
                        interaction2 = 'none'


                    # Overall interaction type

                    if interaction1 == 'positive':
                        if interaction2 == 'positive':
                            overall_interaction = 'mutualism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'none':
                            overall_interaction = 'commensalism'

                    elif interaction1 == 'negative':
                        if interaction2 == 'positive':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'competition'
                        elif interaction2 == 'none':
                            overall_interaction = 'amensalism'


                    elif interaction1 == 'none':
                        if interaction2 == 'positive':
                            overall_interaction = 'commensalism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'amensalism'
                        elif interaction2 == 'none':
                            overall_interaction = 'neutralism'

                    outputting.write(species_1 + ' and ' + species_2 + '\t ' + overall_interaction + '\n')
                    print (species_1 + ' and ' + species_2 + '; ' + overall_interaction + '\n')

            
            print ('Running whole community model.')
            Community_growth = {}

            name1 = models[0].split('/')[-1:][0].split('.RDS')[0]
            name2 = models[1].split('/')[-1:][0].split('.RDS')[0]
            name3 = models[2].split('/')[-1:][0].split('.RDS')[0]
            name4 = models[3].split('/')[-1:][0].split('.RDS')[0]
            name5 = models[4].split('/')[-1:][0].split('.RDS')[0]
            name6 = models[5].split('/')[-1:][0].split('.RDS')[0]
            try:
                subprocess.run(['rm', location_of_call + '/' + 'statistics.csv'])
            except FileNotFoundError:
                lol = 0

            #print( list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1], name1,name2,name3,name4)
            subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size6/Combined_Growth.R', location_of_call + '/' , models[0], models[1],models[2],models[3],models[4],models[5], name1,name2,name3,name4,name5,name6])
            # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4

            ln_num = 0
            growth_info = {}

            Community_growth[models[0].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[1].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[2].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[3].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[4].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[5].split('/')[-1:][0].split('.RDS')[0]] = []

            for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                ln_num +=1
                if ln_num >1:
                    print (",".join(line.split(',')[:-1][1:]))
                    species = line.split(',')[1].replace('"','')
                    growth = int(line.replace('\n','').split(',')[3])
                    Community_growth[species].append(growth) 

            #subprocess.run(['rm ' + location_of_call + '/' + 'statistics.csv'])

            with open(community_outputfile, 'w') as outputting:            
                outputting.write(community_title)

                species_1 = models[0].split('/')[-1:][0].split('.RDS')[0]
                species_2 = models[1].split('/')[-1:][0].split('.RDS')[0]
                species_3 = models[2].split('/')[-1:][0].split('.RDS')[0]
                species_4 = models[3].split('/')[-1:][0].split('.RDS')[0]
                species_5 = models[4].split('/')[-1:][0].split('.RDS')[0]
                species_6 = models[5].split('/')[-1:][0].split('.RDS')[0]


                for species in [species_1,species_2,species_3,species_4,species_5,species_6]:
                    single1 = individual_growths[species][-1:][0]
                    paired1 = Community_growth[species][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    outputting.write(species + '\t' + str(perc1) + '%\n')


            ## Store the data

            all_data_stored[tuple(models)] = [individual_growths,Paired_growths,Community_growth]

            #print ('Completed consortia; ', str(len(all_data_stored.keys())))


        # # Seven species

        # In[ ]:


        if consortia_size_wanted == 7:
            for listi,individual in enumerate(models):
                tmp_models = list(models)
                tmp_models.pop(listi)

                name1 = models[listi].split('/')[-1:][0].split('.RDS')[0]
                name2 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]

                print ('Currently growing; ', name1)


                # !!!!!!!!!!!!! Here, the species being studied is given first, then the other ones, then all the names
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size7/Single_Growth.R', location_of_call + '/' , models[listi], tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5], name1,name2,name3,name4,name5,name6,name7])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                growth_curve = []
                ln_num = 0
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1]
                        growth_curve.append(int(line.replace('\n','').split(',')[3]))
                individual_growths[name1] = growth_curve

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv '])

            print ('Individual growth curves generated for ;' + str(len(individual_growths)))

                
            all_consortia_combinations = list(itertools.combinations(models, 2))

            Paired_growths = {}

            for comb in all_consortia_combinations:
                tmp_models = list(models[:])
                tmp_models.remove(list(comb)[0])
                tmp_models.remove(list(comb)[1])

                name1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                name2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                print ('Strains being grown together;')
                print( list(comb)[0].replace('"','').split('/')[-1:][0].split('.RDS')[0], list(comb)[1].replace('"','').split('/')[-1:][0].split('.RDS')[0], '; names=(', name1,',',name2,')') # Only highlighting those species being paired
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size7/Paired_Growth.R', location_of_call + '/' , list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4], name1,name2,name3,name4,name5,name6,name7])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                ln_num = 0
                growth_info = {}
                growth_info[list(comb)[0].split('/')[-1:][0].split('.RDS')[0]] = []
                growth_info[list(comb)[1].split('/')[-1:][0].split('.RDS')[0]] = []
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1].replace('"','').split('/')[-1:][0].split('.RDS')[0]
                        growth = int(line.replace('\n','').split(',')[3])
                        growth_info[species].append(growth) 
                Paired_growths[comb] = growth_info

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv '])

            ### Identify paired interactions
            with open(paired_outputfile, 'w') as outputting:
                outputting.write(paired_title)
                for comb, data in Paired_growths.items():
                    print ('Calculating interaction between; ')
                    #print (comb)
                    species_1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                    species_2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                    print (species_1)
                    print (species_2)

                    # Calculate interaction type for species1

                    single1 = individual_growths[species_1][-1:][0]
                    paired1 = data[species_1][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    if perc1 >10.0:
                        interaction1 = 'positive'
                    elif perc1 < -10.0:
                        interaction1 = 'negative'
                    else:
                        interaction1 = 'none'


                    # Calculate interaction type for species2

                    single2 = individual_growths[species_2][-1:][0]
                    paired2 = data[species_2][-1:][0]

                    perc2 = ((single2/paired2)*100)-100
                    print (perc1, perc2)

                    if perc2 >10.0:
                        interaction2 = 'positive'
                    elif perc2 < -10.0:
                        interaction2 = 'negative'
                    else:
                        interaction2 = 'none'


                    # Overall interaction type

                    if interaction1 == 'positive':
                        if interaction2 == 'positive':
                            overall_interaction = 'mutualism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'none':
                            overall_interaction = 'commensalism'

                    elif interaction1 == 'negative':
                        if interaction2 == 'positive':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'competition'
                        elif interaction2 == 'none':
                            overall_interaction = 'amensalism'


                    elif interaction1 == 'none':
                        if interaction2 == 'positive':
                            overall_interaction = 'commensalism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'amensalism'
                        elif interaction2 == 'none':
                            overall_interaction = 'neutralism'

                    outputting.write(species_1 + ' and ' + species_2 + '\t' + overall_interaction + '\n')
                    print (species_1 + ' and ' + species_2 + '; ' + overall_interaction + '\n')


            
            print ('Running whole community model.')
            Community_growth = {}

            name1 = models[0].split('/')[-1:][0].split('.RDS')[0]
            name2 = models[1].split('/')[-1:][0].split('.RDS')[0]
            name3 = models[2].split('/')[-1:][0].split('.RDS')[0]
            name4 = models[3].split('/')[-1:][0].split('.RDS')[0]
            name5 = models[4].split('/')[-1:][0].split('.RDS')[0]
            name6 = models[5].split('/')[-1:][0].split('.RDS')[0]
            name7 = models[6].split('/')[-1:][0].split('.RDS')[0]
            try:
                subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])
            except FileNotFoundError:
                lol = 0

            #print( list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1], name1,name2,name3,name4)
            subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size7/Combined_Growth.R', location_of_call + '/' , models[0], models[1],models[2],models[3],models[4],models[5],models[6], name1,name2,name3,name4,name5,name6,name7])
            # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4

            ln_num = 0
            growth_info = {}

            Community_growth[models[0].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[1].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[2].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[3].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[4].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[5].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[6].split('/')[-1:][0].split('.RDS')[0]] = []

            for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                ln_num +=1
                if ln_num >1:
                    print (",".join(line.split(',')[:-1][1:]))
                    species = line.split(',')[1].replace('"','')
                    growth = int(line.replace('\n','').split(',')[3])
                    Community_growth[species].append(growth) 

            #subprocess.run(['rm ' + location_of_call + '/' + 'statistics.csv'])
            with open(community_outputfile, 'w') as outputting:
                outputting.write(community_title)

                species_1 = models[0].split('/')[-1:][0].split('.RDS')[0]
                species_2 = models[1].split('/')[-1:][0].split('.RDS')[0]
                species_3 = models[2].split('/')[-1:][0].split('.RDS')[0]
                species_4 = models[3].split('/')[-1:][0].split('.RDS')[0]
                species_5 = models[4].split('/')[-1:][0].split('.RDS')[0]
                species_6 = models[5].split('/')[-1:][0].split('.RDS')[0]
                species_7 = models[6].split('/')[-1:][0].split('.RDS')[0]


                for species in [species_1,species_2,species_3,species_4,species_5,species_6,species_7]:
                    single1 = individual_growths[species][-1:][0]
                    paired1 = Community_growth[species][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    outputting.write(species + '\t' + str(perc1) + '%\n')

            ## Store the data

            all_data_stored[tuple(models)] = [individual_growths,Paired_growths,Community_growth]

            #print ('Completed consortia; ', str(len(all_data_stored.keys())))


        # # Eight species

        # In[ ]:


        if consortia_size_wanted == 8:
            for listi,individual in enumerate(models):
                tmp_models = list(models)
                tmp_models.pop(listi)

                name1 = models[listi].split('/')[-1:][0].split('.RDS')[0]
                name2 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[6].split('/')[-1:][0].split('.RDS')[0]

                print ('Currently growing; ', name1)


                # !!!!!!!!!!!!! Here, the species being studied is given first, then the other ones, then all the names
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size8/Single_Growth.R', location_of_call + '/' , models[listi], tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5],tmp_models[6], name1,name2,name3,name4,name5,name6,name7,name8])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                growth_curve = []
                ln_num = 0
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1]
                        growth_curve.append(int(line.replace('\n','').split(',')[3]))
                individual_growths[name1] = growth_curve

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv '])

            print ('Individual growth curves generated for ;' + str(len(individual_growths)))

                
            all_consortia_combinations = list(itertools.combinations(models, 2))

            Paired_growths = {}

            for comb in all_consortia_combinations:
                tmp_models = list(models[:])
                tmp_models.remove(list(comb)[0])
                tmp_models.remove(list(comb)[1])

                name1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                name2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                print ('Strains being grown together;')
                print( list(comb)[0].replace('"','').split('/')[-1:][0].split('.RDS')[0], list(comb)[1].replace('"','').split('/')[-1:][0].split('.RDS')[0], '; names=(', name1,',',name2,')') # Only highlighting those species being paired
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size8/Paired_Growth.R', location_of_call + '/' , list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5], name1,name2,name3,name4,name5,name6,name7,name8])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                ln_num = 0
                growth_info = {}
                growth_info[list(comb)[0].split('/')[-1:][0].split('.RDS')[0]] = []
                growth_info[list(comb)[1].split('/')[-1:][0].split('.RDS')[0]] = []
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1].replace('"','').split('/')[-1:][0].split('.RDS')[0]
                        growth = int(line.replace('\n','').split(',')[3])
                        growth_info[species].append(growth) 
                Paired_growths[comb] = growth_info

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv '])

            ### Identify paired interactions
            with open(paired_outputfile, 'w') as outputting:
                outputting.write(paired_title)
                for comb, data in Paired_growths.items():
                    print ('Calculating interaction between; ')
                    #print (comb)
                    species_1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                    species_2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                    print (species_1)
                    print (species_2)

                    # Calculate interaction type for species1

                    single1 = individual_growths[species_1][-1:][0]
                    paired1 = data[species_1][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    if perc1 >10.0:
                        interaction1 = 'positive'
                    elif perc1 < -10.0:
                        interaction1 = 'negative'
                    else:
                        interaction1 = 'none'


                    # Calculate interaction type for species2

                    single2 = individual_growths[species_2][-1:][0]
                    paired2 = data[species_2][-1:][0]

                    perc2 = ((single2/paired2)*100)-100
                    print (perc1, perc2)

                    if perc2 >10.0:
                        interaction2 = 'positive'
                    elif perc2 < -10.0:
                        interaction2 = 'negative'
                    else:
                        interaction2 = 'none'


                    # Overall interaction type

                    if interaction1 == 'positive':
                        if interaction2 == 'positive':
                            overall_interaction = 'mutualism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'none':
                            overall_interaction = 'commensalism'

                    elif interaction1 == 'negative':
                        if interaction2 == 'positive':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'competition'
                        elif interaction2 == 'none':
                            overall_interaction = 'amensalism'


                    elif interaction1 == 'none':
                        if interaction2 == 'positive':
                            overall_interaction = 'commensalism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'amensalism'
                        elif interaction2 == 'none':
                            overall_interaction = 'neutralism'

                    outputting.write(species_1 + ' and ' + species_2 + '\t' + overall_interaction + '\n')
                    print (species_1 + ' and ' + species_2 + '; ' + overall_interaction + '\n')


            
            print ('Running whole community model.')
            Community_growth = {}

            name1 = models[0].split('/')[-1:][0].split('.RDS')[0]
            name2 = models[1].split('/')[-1:][0].split('.RDS')[0]
            name3 = models[2].split('/')[-1:][0].split('.RDS')[0]
            name4 = models[3].split('/')[-1:][0].split('.RDS')[0]
            name5 = models[4].split('/')[-1:][0].split('.RDS')[0]
            name6 = models[5].split('/')[-1:][0].split('.RDS')[0]
            name7 = models[6].split('/')[-1:][0].split('.RDS')[0]
            name8 = models[7].split('/')[-1:][0].split('.RDS')[0]
            try:
                subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])
            except FileNotFoundError:
                lol = 0

            #print( list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1], name1,name2,name3,name4)
            subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size8/Combined_Growth.R', location_of_call + '/' , models[0], models[1],models[2],models[3],models[4],models[5],models[6],models[7], name1,name2,name3,name4,name5,name6,name7,name8])
            # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4

            ln_num = 0
            growth_info = {}

            Community_growth[models[0].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[1].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[2].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[3].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[4].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[5].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[6].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[7].split('/')[-1:][0].split('.RDS')[0]] = []

            for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                ln_num +=1
                if ln_num >1:
                    print (",".join(line.split(',')[:-1][1:]))
                    species = line.split(',')[1].replace('"','')
                    growth = int(line.replace('\n','').split(',')[3])
                    Community_growth[species].append(growth) 

            #subprocess.run(['rm ' + location_of_call + '/' + 'statistics.csv'])
            with open(community_outputfile, 'w') as outputting:
                outputting.write(community_title)

                species_1 = models[0].split('/')[-1:][0].split('.RDS')[0]
                species_2 = models[1].split('/')[-1:][0].split('.RDS')[0]
                species_3 = models[2].split('/')[-1:][0].split('.RDS')[0]
                species_4 = models[3].split('/')[-1:][0].split('.RDS')[0]
                species_5 = models[4].split('/')[-1:][0].split('.RDS')[0]
                species_6 = models[5].split('/')[-1:][0].split('.RDS')[0]
                species_7 = models[6].split('/')[-1:][0].split('.RDS')[0]
                species_8 = models[7].split('/')[-1:][0].split('.RDS')[0]


                for species in [species_1,species_2,species_3,species_4,species_5,species_6,species_7,species_8]:
                    single1 = individual_growths[species][-1:][0]
                    paired1 = Community_growth[species][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    outputting.write(species + '\t' + str(perc1) + '%\n')

            ## Store the data

            all_data_stored[tuple(models)] = [individual_growths,Paired_growths,Community_growth]

            #print ('Completed consortia; ', str(len(all_data_stored.keys())))


        # # Nine species

        # In[ ]:


        if consortia_size_wanted == 9:
            for listi,individual in enumerate(models):
                tmp_models = list(models)
                tmp_models.pop(listi)

                name1 = models[listi].split('/')[-1:][0].split('.RDS')[0]
                name2 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[6].split('/')[-1:][0].split('.RDS')[0]
                name9 = tmp_models[7].split('/')[-1:][0].split('.RDS')[0]

                print ('Currently growing; ', name1)


                # !!!!!!!!!!!!! Here, the species being studied is given first, then the other ones, then all the names
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size9/Single_Growth.R', location_of_call + '/' , models[listi], tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5],tmp_models[6],tmp_models[7], name1,name2,name3,name4,name5,name6,name7,name8,name9])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                growth_curve = []
                ln_num = 0
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1]
                        growth_curve.append(int(line.replace('\n','').split(',')[3]))
                individual_growths[name1] = growth_curve

            subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])

            print ('Individual growth curves generated for ;' + str(len(individual_growths)))

                
            all_consortia_combinations = list(itertools.combinations(models, 2))

            Paired_growths = {}

            for comb in all_consortia_combinations:
                tmp_models = list(models[:])
                tmp_models.remove(list(comb)[0])
                tmp_models.remove(list(comb)[1])

                name1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                name2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                name9 = tmp_models[6].split('/')[-1:][0].split('.RDS')[0]
                print ('Strains being grown together;')
                print( list(comb)[0].replace('"','').split('/')[-1:][0].split('.RDS')[0], list(comb)[1].replace('"','').split('/')[-1:][0].split('.RDS')[0], '; names=(', name1,',',name2,')') # Only highlighting those species being paired
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size9/Paired_Growth.R', location_of_call + '/' , list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5],tmp_models[6], name1,name2,name3,name4,name5,name6,name7,name8,name9])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                ln_num = 0
                growth_info = {}
                growth_info[list(comb)[0].split('/')[-1:][0].split('.RDS')[0]] = []
                growth_info[list(comb)[1].split('/')[-1:][0].split('.RDS')[0]] = []
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1].replace('"','').split('/')[-1:][0].split('.RDS')[0]
                        growth = int(line.replace('\n','').split(',')[3])
                        growth_info[species].append(growth) 
                Paired_growths[comb] = growth_info

            subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])

            ### Identify paired interactions
            with open(paired_outputfile, 'w') as outputting:
                outputting.write(paired_title)
                for comb, data in Paired_growths.items():
                    print ('Calculating interaction between; ')
                    #print (comb)
                    species_1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                    species_2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                    print (species_1)
                    print (species_2)

                    # Calculate interaction type for species1

                    single1 = individual_growths[species_1][-1:][0]
                    paired1 = data[species_1][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    if perc1 >10.0:
                        interaction1 = 'positive'
                    elif perc1 < -10.0:
                        interaction1 = 'negative'
                    else:
                        interaction1 = 'none'


                    # Calculate interaction type for species2

                    single2 = individual_growths[species_2][-1:][0]
                    paired2 = data[species_2][-1:][0]

                    perc2 = ((single2/paired2)*100)-100
                    print (perc1, perc2)

                    if perc2 >10.0:
                        interaction2 = 'positive'
                    elif perc2 < -10.0:
                        interaction2 = 'negative'
                    else:
                        interaction2 = 'none'


                    # Overall interaction type

                    if interaction1 == 'positive':
                        if interaction2 == 'positive':
                            overall_interaction = 'mutualism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'none':
                            overall_interaction = 'commensalism'

                    elif interaction1 == 'negative':
                        if interaction2 == 'positive':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'competition'
                        elif interaction2 == 'none':
                            overall_interaction = 'amensalism'


                    elif interaction1 == 'none':
                        if interaction2 == 'positive':
                            overall_interaction = 'commensalism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'amensalism'
                        elif interaction2 == 'none':
                            overall_interaction = 'neutralism'

                    outputting.write(species_1 + ' and ' + species_2 + '\t' + overall_interaction + '\n')
                    print (species_1 + ' and ' + species_2 + '; ' + overall_interaction + '\n')

            
            print ('Running whole community model.')
            Community_growth = {}

            name1 = models[0].split('/')[-1:][0].split('.RDS')[0]
            name2 = models[1].split('/')[-1:][0].split('.RDS')[0]
            name3 = models[2].split('/')[-1:][0].split('.RDS')[0]
            name4 = models[3].split('/')[-1:][0].split('.RDS')[0]
            name5 = models[4].split('/')[-1:][0].split('.RDS')[0]
            name6 = models[5].split('/')[-1:][0].split('.RDS')[0]
            name7 = models[6].split('/')[-1:][0].split('.RDS')[0]
            name8 = models[7].split('/')[-1:][0].split('.RDS')[0]
            name9 = models[8].split('/')[-1:][0].split('.RDS')[0]
            try:
                subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])
            except FileNotFoundError:
                lol = 0

            #print( list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1], name1,name2,name3,name4)
            subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size9/Combined_Growth.R', location_of_call + '/' , models[0], models[1],models[2],models[3],models[4],models[5],models[6],models[7],models[8], name1,name2,name3,name4,name5,name6,name7,name8,name9])
            # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4

            ln_num = 0
            growth_info = {}

            Community_growth[models[0].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[1].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[2].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[3].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[4].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[5].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[6].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[7].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[8].split('/')[-1:][0].split('.RDS')[0]] = []

            for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                ln_num +=1
                if ln_num >1:
                    print (",".join(line.split(',')[:-1][1:]))
                    species = line.split(',')[1].replace('"','')
                    growth = int(line.replace('\n','').split(',')[3])
                    Community_growth[species].append(growth) 

            #subprocess.run(['rm ' + location_of_call + '/' + 'statistics.csv'])
            with open(community_outputfile, 'w') as outputting:
                outputting.write(community_title)

                species_1 = models[0].split('/')[-1:][0].split('.RDS')[0]
                species_2 = models[1].split('/')[-1:][0].split('.RDS')[0]
                species_3 = models[2].split('/')[-1:][0].split('.RDS')[0]
                species_4 = models[3].split('/')[-1:][0].split('.RDS')[0]
                species_5 = models[4].split('/')[-1:][0].split('.RDS')[0]
                species_6 = models[5].split('/')[-1:][0].split('.RDS')[0]
                species_7 = models[6].split('/')[-1:][0].split('.RDS')[0]
                species_8 = models[7].split('/')[-1:][0].split('.RDS')[0]
                species_9 = models[8].split('/')[-1:][0].split('.RDS')[0]


                for species in [species_1,species_2,species_3,species_4,species_5,species_6,species_7,species_8,species_9]:
                    single1 = individual_growths[species][-1:][0]
                    paired1 = Community_growth[species][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    outputting.write(species + '\t' + str(perc1) + '%\n')

            outputting.close()


            ## Store the data

            all_data_stored[tuple(models)] = [individual_growths,Paired_growths,Community_growth]

            #print ('Completed consortia; ', str(len(all_data_stored.keys())))


        # In[105]:


        # Ten species


        # In[66]:


        if consortia_size_wanted == 10:
            for listi,individual in enumerate(models):
                tmp_models = list(models)
                tmp_models.pop(listi)

                name1 = models[listi].split('/')[-1:][0].split('.RDS')[0]
                name2 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[6].split('/')[-1:][0].split('.RDS')[0]
                name9 = tmp_models[7].split('/')[-1:][0].split('.RDS')[0]
                name10 = tmp_models[8].split('/')[-1:][0].split('.RDS')[0]

                print ('Currently growing; ', name1)

                # !!!!!!!!!!!!! Here, the species being studied is given first, then the other ones, then all the names
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size10/Single_Growth.R', location_of_call + '/' , models[listi], tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5],tmp_models[6],tmp_models[7],tmp_models[8], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                growth_curve = []
                ln_num = 0

                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1]
                        growth_curve.append(int(line.replace('\n','').split(',')[3]))
                individual_growths[name1] = growth_curve

            subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])

            print ('Individual growth curves generated for ;' + str(len(individual_growths)))

                
            all_consortia_combinations = list(itertools.combinations(models, 2))

            Paired_growths = {}

            for comb in all_consortia_combinations:
                tmp_models = list(models[:])
                tmp_models.remove(list(comb)[0])
                tmp_models.remove(list(comb)[1])

                name1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                name2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                name9 = tmp_models[6].split('/')[-1:][0].split('.RDS')[0]
                name10 = tmp_models[7].split('/')[-1:][0].split('.RDS')[0]
                print ('Strains being grown together;')
                print( list(comb)[0].replace('"','').split('/')[-1:][0].split('.RDS')[0], list(comb)[1].replace('"','').split('/')[-1:][0].split('.RDS')[0], '; names=(', name1,',',name2,')') # Only highlighting those species being paired
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size10/Paired_Growth.R', location_of_call + '/' , list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5],tmp_models[6],tmp_models[7], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                ln_num = 0
                growth_info = {}
                growth_info[list(comb)[0].split('/')[-1:][0].split('.RDS')[0]] = []
                growth_info[list(comb)[1].split('/')[-1:][0].split('.RDS')[0]] = []
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1].replace('"','').split('/')[-1:][0].split('.RDS')[0]
                        growth = int(line.replace('\n','').split(',')[3])
                        growth_info[species].append(growth) 
                Paired_growths[comb] = growth_info

            subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])

            ### Identify paired interactions
            with open(paired_outputfile, 'w') as outputting:
                outputting.write(paired_title)
                for comb, data in Paired_growths.items():
                    print ('Calculating interaction between; ')
                    #print (comb)
                    species_1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                    species_2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                    print (species_1)
                    print (species_2)

                    # Calculate interaction type for species1

                    single1 = individual_growths[species_1][-1:][0]
                    paired1 = data[species_1][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    if perc1 >10.0:
                        interaction1 = 'positive'
                    elif perc1 < -10.0:
                        interaction1 = 'negative'
                    else:
                        interaction1 = 'none'


                    # Calculate interaction type for species2

                    single2 = individual_growths[species_2][-1:][0]
                    paired2 = data[species_2][-1:][0]

                    perc2 = ((single2/paired2)*100)-100
                    print (perc1, perc2)

                    if perc2 >10.0:
                        interaction2 = 'positive'
                    elif perc2 < -10.0:
                        interaction2 = 'negative'
                    else:
                        interaction2 = 'none'


                    # Overall interaction type

                    if interaction1 == 'positive':
                        if interaction2 == 'positive':
                            overall_interaction = 'mutualism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'none':
                            overall_interaction = 'commensalism'

                    elif interaction1 == 'negative':
                        if interaction2 == 'positive':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'competition'
                        elif interaction2 == 'none':
                            overall_interaction = 'amensalism'


                    elif interaction1 == 'none':
                        if interaction2 == 'positive':
                            overall_interaction = 'commensalism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'amensalism'
                        elif interaction2 == 'none':
                            overall_interaction = 'neutralism'

                    outputting.write(species_1 + ' and ' + species_2 + '\t' + overall_interaction + '\n')
                    print (species_1 + ' and ' + species_2 + '; ' + overall_interaction + '\n')


            
            print ('Running whole community model.')
            Community_growth = {}

            name1 = models[0].split('/')[-1:][0].split('.RDS')[0]
            name2 = models[1].split('/')[-1:][0].split('.RDS')[0]
            name3 = models[2].split('/')[-1:][0].split('.RDS')[0]
            name4 = models[3].split('/')[-1:][0].split('.RDS')[0]
            name5 = models[4].split('/')[-1:][0].split('.RDS')[0]
            name6 = models[5].split('/')[-1:][0].split('.RDS')[0]
            name7 = models[6].split('/')[-1:][0].split('.RDS')[0]
            name8 = models[7].split('/')[-1:][0].split('.RDS')[0]
            name9 = models[8].split('/')[-1:][0].split('.RDS')[0]
            name10 = models[9].split('/')[-1:][0].split('.RDS')[0]
            try:
                subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])
            except FileNotFoundError:
                lol = 0

            #print( list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1], name1,name2,name3,name4)
            subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size10/Combined_Growth.R', location_of_call + '/' , models[0], models[1],models[2],models[3],models[4],models[5],models[6],models[7],models[8],models[9], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10])
            # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4

            ln_num = 0
            growth_info = {}

            Community_growth[models[0].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[1].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[2].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[3].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[4].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[5].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[6].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[7].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[8].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[9].split('/')[-1:][0].split('.RDS')[0]] = []

            for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                ln_num +=1
                if ln_num >1:
                    print (",".join(line.split(',')[:-1][1:]))
                    species = line.split(',')[1].replace('"','')
                    growth = int(line.replace('\n','').split(',')[3])
                    Community_growth[species].append(growth) 

            #subprocess.run(['rm ' + location_of_call + '/' + 'statistics.csv'])
            with open(community_outputfile, 'w') as outputting:
                outputting.write(community_title)

                species_1 = models[0].split('/')[-1:][0].split('.RDS')[0]
                species_2 = models[1].split('/')[-1:][0].split('.RDS')[0]
                species_3 = models[2].split('/')[-1:][0].split('.RDS')[0]
                species_4 = models[3].split('/')[-1:][0].split('.RDS')[0]
                species_5 = models[4].split('/')[-1:][0].split('.RDS')[0]
                species_6 = models[5].split('/')[-1:][0].split('.RDS')[0]
                species_7 = models[6].split('/')[-1:][0].split('.RDS')[0]
                species_8 = models[7].split('/')[-1:][0].split('.RDS')[0]
                species_9 = models[8].split('/')[-1:][0].split('.RDS')[0]
                species_10 = models[9].split('/')[-1:][0].split('.RDS')[0]


                for species in [species_1,species_2,species_3,species_4,species_5,species_6,species_7,species_8,species_9,species_10]:
                    single1 = individual_growths[species][-1:][0]
                    paired1 = Community_growth[species][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    outputting.write(species + '\t' + str(perc1) + '%\n')



            ## Store the data

            all_data_stored[tuple(models)] = [individual_growths,Paired_growths,Community_growth]

            #print ('Completed consortia; ', str(len(all_data_stored.keys())))


        # In[106]:


        # Eleven species


        # In[ ]:


        if consortia_size_wanted == 11:
            for listi,individual in enumerate(models):
                tmp_models = list(models)
                tmp_models.pop(listi)

                name1 = models[listi].split('/')[-1:][0].split('.RDS')[0]
                name2 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[6].split('/')[-1:][0].split('.RDS')[0]
                name9 = tmp_models[7].split('/')[-1:][0].split('.RDS')[0]
                name10 = tmp_models[8].split('/')[-1:][0].split('.RDS')[0]
                name11 = tmp_models[9].split('/')[-1:][0].split('.RDS')[0]

                print ('Currently growing; ', name1)


                # !!!!!!!!!!!!! Here, the species being studied is given first, then the other ones, then all the names
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size11/Single_Growth.R', location_of_call + '/' , models[listi], tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5],tmp_models[6],tmp_models[7],tmp_models[8],tmp_models[9], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                growth_curve = []
                ln_num = 0
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1]
                        growth_curve.append(int(line.replace('\n','').split(',')[3]))
                individual_growths[name1] = growth_curve

            subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])

            print ('Individual growth curves generated for ;' + str(len(individual_growths)))

                
            all_consortia_combinations = list(itertools.combinations(models, 2))

            Paired_growths = {}

            for comb in all_consortia_combinations:
                tmp_models = list(models[:])
                tmp_models.remove(list(comb)[0])
                tmp_models.remove(list(comb)[1])

                name1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                name2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                name9 = tmp_models[6].split('/')[-1:][0].split('.RDS')[0]
                name10 = tmp_models[7].split('/')[-1:][0].split('.RDS')[0]
                name11 = tmp_models[8].split('/')[-1:][0].split('.RDS')[0]
                print ('Strains being grown together;')
                print( list(comb)[0].replace('"','').split('/')[-1:][0].split('.RDS')[0], list(comb)[1].replace('"','').split('/')[-1:][0].split('.RDS')[0], '; names=(', name1,',',name2,')') # Only highlighting those species being paired
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size11/Paired_Growth.R', location_of_call + '/' , list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5],tmp_models[6],tmp_models[7],tmp_models[8], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                ln_num = 0
                growth_info = {}
                growth_info[list(comb)[0].split('/')[-1:][0].split('.RDS')[0]] = []
                growth_info[list(comb)[1].split('/')[-1:][0].split('.RDS')[0]] = []
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1].replace('"','').split('/')[-1:][0].split('.RDS')[0]
                        growth = int(line.replace('\n','').split(',')[3])
                        growth_info[species].append(growth) 
                Paired_growths[comb] = growth_info

            subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])

            ### Identify paired interactions
            with open(paired_outputfile, 'w') as outputting:
                outputting.write(paired_title)
                for comb, data in Paired_growths.items():
                    print ('Calculating interaction between; ')
                    #print (comb)
                    species_1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                    species_2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                    print (species_1)
                    print (species_2)

                    # Calculate interaction type for species1

                    single1 = individual_growths[species_1][-1:][0]
                    paired1 = data[species_1][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    if perc1 >10.0:
                        interaction1 = 'positive'
                    elif perc1 < -10.0:
                        interaction1 = 'negative'
                    else:
                        interaction1 = 'none'


                    # Calculate interaction type for species2

                    single2 = individual_growths[species_2][-1:][0]
                    paired2 = data[species_2][-1:][0]

                    perc2 = ((single2/paired2)*100)-100
                    print (perc1, perc2)

                    if perc2 >10.0:
                        interaction2 = 'positive'
                    elif perc2 < -10.0:
                        interaction2 = 'negative'
                    else:
                        interaction2 = 'none'


                    # Overall interaction type

                    if interaction1 == 'positive':
                        if interaction2 == 'positive':
                            overall_interaction = 'mutualism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'none':
                            overall_interaction = 'commensalism'

                    elif interaction1 == 'negative':
                        if interaction2 == 'positive':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'competition'
                        elif interaction2 == 'none':
                            overall_interaction = 'amensalism'


                    elif interaction1 == 'none':
                        if interaction2 == 'positive':
                            overall_interaction = 'commensalism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'amensalism'
                        elif interaction2 == 'none':
                            overall_interaction = 'neutralism'

                    outputting.write(species_1 + ' and ' + species_2 + '\t' + overall_interaction + '\n')
                    print (species_1 + ' and ' + species_2 + '; ' + overall_interaction + '\n')

            
            print ('Running whole community model.')
            Community_growth = {}

            name1 = models[0].split('/')[-1:][0].split('.RDS')[0]
            name2 = models[1].split('/')[-1:][0].split('.RDS')[0]
            name3 = models[2].split('/')[-1:][0].split('.RDS')[0]
            name4 = models[3].split('/')[-1:][0].split('.RDS')[0]
            name5 = models[4].split('/')[-1:][0].split('.RDS')[0]
            name6 = models[5].split('/')[-1:][0].split('.RDS')[0]
            name7 = models[6].split('/')[-1:][0].split('.RDS')[0]
            name8 = models[7].split('/')[-1:][0].split('.RDS')[0]
            name9 = models[8].split('/')[-1:][0].split('.RDS')[0]
            name10 = models[9].split('/')[-1:][0].split('.RDS')[0]
            name11 = models[10].split('/')[-1:][0].split('.RDS')[0]
            try:
                subprocess.run(['rm', location_of_call + '/' + 'statistics.csv'])
            except FileNotFoundError:
                lol = 0

            #print( list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1], name1,name2,name3,name4)
            subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size11/Combined_Growth.R', location_of_call + '/' , models[0], models[1],models[2],models[3],models[4],models[5],models[6],models[7],models[8],models[9],models[10], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11])
            # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4

            ln_num = 0
            growth_info = {}

            Community_growth[models[0].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[1].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[2].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[3].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[4].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[5].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[6].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[7].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[8].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[9].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[10].split('/')[-1:][0].split('.RDS')[0]] = []

            for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                ln_num +=1
                if ln_num >1:
                    print (",".join(line.split(',')[:-1][1:]))
                    species = line.split(',')[1].replace('"','')
                    growth = int(line.replace('\n','').split(',')[3])
                    Community_growth[species].append(growth) 

            #subprocess.run(['rm ' + location_of_call + '/' + 'statistics.csv'])
            with open(community_outputfile, 'w') as outputting:
                outputting.write(community_title)

                species_1 = models[0].split('/')[-1:][0].split('.RDS')[0]
                species_2 = models[1].split('/')[-1:][0].split('.RDS')[0]
                species_3 = models[2].split('/')[-1:][0].split('.RDS')[0]
                species_4 = models[3].split('/')[-1:][0].split('.RDS')[0]
                species_5 = models[4].split('/')[-1:][0].split('.RDS')[0]
                species_6 = models[5].split('/')[-1:][0].split('.RDS')[0]
                species_7 = models[6].split('/')[-1:][0].split('.RDS')[0]
                species_8 = models[7].split('/')[-1:][0].split('.RDS')[0]
                species_9 = models[8].split('/')[-1:][0].split('.RDS')[0]
                species_10 = models[9].split('/')[-1:][0].split('.RDS')[0]
                species_11 = models[10].split('/')[-1:][0].split('.RDS')[0]


                for species in [species_1,species_2,species_3,species_4,species_5,species_6,species_7,species_8,species_9,species_10,species_11]:
                    single1 = individual_growths[species][-1:][0]
                    paired1 = Community_growth[species][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    outputting.write(species + '\t' + str(perc1) + '%\n')



            ## Store the data

            all_data_stored[tuple(models)] = [individual_growths,Paired_growths,Community_growth]

            #print ('Completed consortia; ', str(len(all_data_stored.keys())))


        # In[107]:


        # Twelve species


        # In[ ]:


        if consortia_size_wanted == 12:
            for listi,individual in enumerate(models):
                tmp_models = list(models)
                tmp_models.pop(listi)

                name1 = models[listi].split('/')[-1:][0].split('.RDS')[0]
                name2 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[6].split('/')[-1:][0].split('.RDS')[0]
                name9 = tmp_models[7].split('/')[-1:][0].split('.RDS')[0]
                name10 = tmp_models[8].split('/')[-1:][0].split('.RDS')[0]
                name11 = tmp_models[9].split('/')[-1:][0].split('.RDS')[0]
                name12 = tmp_models[10].split('/')[-1:][0].split('.RDS')[0]

                print ('Currently growing; ', name1)


                # !!!!!!!!!!!!! Here, the species being studied is given first, then the other ones, then all the names
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size11/Single_Growth.R', location_of_call + '/' , models[listi], tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5],tmp_models[6],tmp_models[7],tmp_models[8],tmp_models[9],tmp_models[10], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11,name12])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                growth_curve = []
                ln_num = 0
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1]
                        growth_curve.append(int(line.replace('\n','').split(',')[3]))
                individual_growths[name1] = growth_curve

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv'])

            print ('Individual growth curves generated for ;' + str(len(individual_growths)))

                
            all_consortia_combinations = list(itertools.combinations(models, 2))

            Paired_growths = {}

            for comb in all_consortia_combinations:
                tmp_models = list(models[:])
                tmp_models.remove(list(comb)[0])
                tmp_models.remove(list(comb)[1])

                name1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                name2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                name9 = tmp_models[6].split('/')[-1:][0].split('.RDS')[0]
                name10 = tmp_models[7].split('/')[-1:][0].split('.RDS')[0]
                name11 = tmp_models[8].split('/')[-1:][0].split('.RDS')[0]
                name12 = tmp_models[9].split('/')[-1:][0].split('.RDS')[0]
                print ('Strains being grown together;')
                print( list(comb)[0].replace('"','').split('/')[-1:][0].split('.RDS')[0], list(comb)[1].replace('"','').split('/')[-1:][0].split('.RDS')[0], '; names=(', name1,',',name2,')') # Only highlighting those species being paired
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size12/Paired_Growth.R', location_of_call + '/' , list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5],tmp_models[6],tmp_models[7],tmp_models[8],tmp_models[9], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11,name12])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                ln_num = 0
                growth_info = {}
                growth_info[list(comb)[0].split('/')[-1:][0].split('.RDS')[0]] = []
                growth_info[list(comb)[1].split('/')[-1:][0].split('.RDS')[0]] = []
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1].replace('"','').split('/')[-1:][0].split('.RDS')[0]
                        growth = int(line.replace('\n','').split(',')[3])
                        growth_info[species].append(growth) 
                Paired_growths[comb] = growth_info

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv'])

            ### Identify paired interactions
            with open(paired_outputfile, 'w') as outputting:
                outputting.write(paired_title)
                for comb, data in Paired_growths.items():
                    print ('Calculating interaction between; ')
                    #print (comb)
                    species_1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                    species_2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                    print (species_1)
                    print (species_2)

                    # Calculate interaction type for species1

                    single1 = individual_growths[species_1][-1:][0]
                    paired1 = data[species_1][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    if perc1 >10.0:
                        interaction1 = 'positive'
                    elif perc1 < -10.0:
                        interaction1 = 'negative'
                    else:
                        interaction1 = 'none'


                    # Calculate interaction type for species2

                    single2 = individual_growths[species_2][-1:][0]
                    paired2 = data[species_2][-1:][0]

                    perc2 = ((single2/paired2)*100)-100
                    print (perc1, perc2)

                    if perc2 >10.0:
                        interaction2 = 'positive'
                    elif perc2 < -10.0:
                        interaction2 = 'negative'
                    else:
                        interaction2 = 'none'


                    # Overall interaction type

                    if interaction1 == 'positive':
                        if interaction2 == 'positive':
                            overall_interaction = 'mutualism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'none':
                            overall_interaction = 'commensalism'

                    elif interaction1 == 'negative':
                        if interaction2 == 'positive':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'competition'
                        elif interaction2 == 'none':
                            overall_interaction = 'amensalism'


                    elif interaction1 == 'none':
                        if interaction2 == 'positive':
                            overall_interaction = 'commensalism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'amensalism'
                        elif interaction2 == 'none':
                            overall_interaction = 'neutralism'

                    outputting.write(species_1 + ' and ' + species_2 + '\t' + overall_interaction + '\n')
                    print (species_1 + ' and ' + species_2 + '; ' + overall_interaction + '\n')

            
            print ('Running whole community model.')
            Community_growth = {}

            name1 = models[0].split('/')[-1:][0].split('.RDS')[0]
            name2 = models[1].split('/')[-1:][0].split('.RDS')[0]
            name3 = models[2].split('/')[-1:][0].split('.RDS')[0]
            name4 = models[3].split('/')[-1:][0].split('.RDS')[0]
            name5 = models[4].split('/')[-1:][0].split('.RDS')[0]
            name6 = models[5].split('/')[-1:][0].split('.RDS')[0]
            name7 = models[6].split('/')[-1:][0].split('.RDS')[0]
            name8 = models[7].split('/')[-1:][0].split('.RDS')[0]
            name9 = models[8].split('/')[-1:][0].split('.RDS')[0]
            name10 = models[9].split('/')[-1:][0].split('.RDS')[0]
            name11 = models[10].split('/')[-1:][0].split('.RDS')[0]
            name12 = models[11].split('/')[-1:][0].split('.RDS')[0]
            try:
                subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])
            except FileNotFoundError:
                lol = 0

            #print( list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1], name1,name2,name3,name4)
            subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size12/Combined_Growth.R', location_of_call + '/' , models[0], models[1],models[2],models[3],models[4],models[5],models[6],models[7],models[8],models[9],models[10],models[11], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11,name12])
            # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4

            ln_num = 0
            growth_info = {}

            Community_growth[models[0].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[1].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[2].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[3].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[4].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[5].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[6].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[7].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[8].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[9].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[10].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[11].split('/')[-1:][0].split('.RDS')[0]] = []

            for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                ln_num +=1
                if ln_num >1:
                    print (",".join(line.split(',')[:-1][1:]))
                    species = line.split(',')[1].replace('"','')
                    growth = int(line.replace('\n','').split(',')[3])
                    Community_growth[species].append(growth) 

            #subprocess.run(['rm ' + location_of_call + '/' + 'statistics.csv'])
            with open(community_outputfile, 'w') as outputting:
                outputting.write(community_title)

                species_1 = models[0].split('/')[-1:][0].split('.RDS')[0]
                species_2 = models[1].split('/')[-1:][0].split('.RDS')[0]
                species_3 = models[2].split('/')[-1:][0].split('.RDS')[0]
                species_4 = models[3].split('/')[-1:][0].split('.RDS')[0]
                species_5 = models[4].split('/')[-1:][0].split('.RDS')[0]
                species_6 = models[5].split('/')[-1:][0].split('.RDS')[0]
                species_7 = models[6].split('/')[-1:][0].split('.RDS')[0]
                species_8 = models[7].split('/')[-1:][0].split('.RDS')[0]
                species_9 = models[8].split('/')[-1:][0].split('.RDS')[0]
                species_10 = models[9].split('/')[-1:][0].split('.RDS')[0]
                species_11 = models[10].split('/')[-1:][0].split('.RDS')[0]
                species_12 = models[11].split('/')[-1:][0].split('.RDS')[0]


                for species in [species_1,species_2,species_3,species_4,species_5,species_6,species_7,species_8,species_9,species_10,species_11,species_12]:
                    single1 = individual_growths[species][-1:][0]
                    paired1 = Community_growth[species][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    outputting.write(species + '\t' + str(perc1) + '%\n')



            ## Store the data

            all_data_stored[tuple(models)] = [individual_growths,Paired_growths,Community_growth]

            #print ('Completed consortia; ', str(len(all_data_stored.keys())))


        # In[108]:


        # Thirteen species


        # In[ ]:


        if consortia_size_wanted == 13:
            for listi,individual in enumerate(models):
                tmp_models = list(models)
                tmp_models.pop(listi)

                name1 = models[listi].split('/')[-1:][0].split('.RDS')[0]
                name2 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[6].split('/')[-1:][0].split('.RDS')[0]
                name9 = tmp_models[7].split('/')[-1:][0].split('.RDS')[0]
                name10 = tmp_models[8].split('/')[-1:][0].split('.RDS')[0]
                name11 = tmp_models[9].split('/')[-1:][0].split('.RDS')[0]
                name12 = tmp_models[10].split('/')[-1:][0].split('.RDS')[0]
                name13 = tmp_models[11].split('/')[-1:][0].split('.RDS')[0]

                print ('Currently growing; ', name1)


                # !!!!!!!!!!!!! Here, the species being studied is given first, then the other ones, then all the names
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size13/Single_Growth.R', location_of_call + '/' , models[listi], tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5],tmp_models[6],tmp_models[7],tmp_models[8],tmp_models[9],tmp_models[10],tmp_models[11], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11,name12,name13])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                growth_curve = []
                ln_num = 0
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1]
                        growth_curve.append(int(line.replace('\n','').split(',')[3]))
                individual_growths[name1] = growth_curve

            subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])

            print ('Individual growth curves generated for ;' + str(len(individual_growths)))

                
            all_consortia_combinations = list(itertools.combinations(models, 2))

            Paired_growths = {}

            for comb in all_consortia_combinations:
                tmp_models = list(models[:])
                tmp_models.remove(list(comb)[0])
                tmp_models.remove(list(comb)[1])

                name1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                name2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                name9 = tmp_models[6].split('/')[-1:][0].split('.RDS')[0]
                name10 = tmp_models[7].split('/')[-1:][0].split('.RDS')[0]
                name11 = tmp_models[8].split('/')[-1:][0].split('.RDS')[0]
                name12 = tmp_models[9].split('/')[-1:][0].split('.RDS')[0]
                name13 = tmp_models[10].split('/')[-1:][0].split('.RDS')[0]
                print ('Strains being grown together;')
                print( list(comb)[0].replace('"','').split('/')[-1:][0].split('.RDS')[0], list(comb)[1].replace('"','').split('/')[-1:][0].split('.RDS')[0], '; names=(', name1,',',name2,')') # Only highlighting those species being paired
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size13/Paired_Growth.R', location_of_call + '/' , list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5],tmp_models[6],tmp_models[7],tmp_models[8],tmp_models[9],tmp_models[10], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11,name12,name13])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                ln_num = 0
                growth_info = {}
                growth_info[list(comb)[0].split('/')[-1:][0].split('.RDS')[0]] = []
                growth_info[list(comb)[1].split('/')[-1:][0].split('.RDS')[0]] = []
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1].replace('"','').split('/')[-1:][0].split('.RDS')[0]
                        growth = int(line.replace('\n','').split(',')[3])
                        growth_info[species].append(growth) 
                Paired_growths[comb] = growth_info

            subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])

            ### Identify paired interactions
            with open(paired_outputfile, 'w') as outputting:
                outputting.write(paired_title)
                for comb, data in Paired_growths.items():
                    print ('Calculating interaction between; ')
                    #print (comb)
                    species_1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                    species_2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                    print (species_1)
                    print (species_2)

                    # Calculate interaction type for species1

                    single1 = individual_growths[species_1][-1:][0]
                    paired1 = data[species_1][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    if perc1 >10.0:
                        interaction1 = 'positive'
                    elif perc1 < -10.0:
                        interaction1 = 'negative'
                    else:
                        interaction1 = 'none'


                    # Calculate interaction type for species2

                    single2 = individual_growths[species_2][-1:][0]
                    paired2 = data[species_2][-1:][0]

                    perc2 = ((single2/paired2)*100)-100
                    print (perc1, perc2)

                    if perc2 >10.0:
                        interaction2 = 'positive'
                    elif perc2 < -10.0:
                        interaction2 = 'negative'
                    else:
                        interaction2 = 'none'


                    # Overall interaction type

                    if interaction1 == 'positive':
                        if interaction2 == 'positive':
                            overall_interaction = 'mutualism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'none':
                            overall_interaction = 'commensalism'

                    elif interaction1 == 'negative':
                        if interaction2 == 'positive':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'competition'
                        elif interaction2 == 'none':
                            overall_interaction = 'amensalism'


                    elif interaction1 == 'none':
                        if interaction2 == 'positive':
                            overall_interaction = 'commensalism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'amensalism'
                        elif interaction2 == 'none':
                            overall_interaction = 'neutralism'

                    outputting.write(species_1 + ' and ' + species_2 + '\t' + overall_interaction + '\n')
                    print (species_1 + ' and ' + species_2 + '; ' + overall_interaction + '\n')

            
            print ('Running whole community model.')
            Community_growth = {}

            name1 = models[0].split('/')[-1:][0].split('.RDS')[0]
            name2 = models[1].split('/')[-1:][0].split('.RDS')[0]
            name3 = models[2].split('/')[-1:][0].split('.RDS')[0]
            name4 = models[3].split('/')[-1:][0].split('.RDS')[0]
            name5 = models[4].split('/')[-1:][0].split('.RDS')[0]
            name6 = models[5].split('/')[-1:][0].split('.RDS')[0]
            name7 = models[6].split('/')[-1:][0].split('.RDS')[0]
            name8 = models[7].split('/')[-1:][0].split('.RDS')[0]
            name9 = models[8].split('/')[-1:][0].split('.RDS')[0]
            name10 = models[9].split('/')[-1:][0].split('.RDS')[0]
            name11 = models[10].split('/')[-1:][0].split('.RDS')[0]
            name12 = models[11].split('/')[-1:][0].split('.RDS')[0]
            name13 = models[12].split('/')[-1:][0].split('.RDS')[0]
            try:
                subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])
            except FileNotFoundError:
                lol = 0

            #print( list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1], name1,name2,name3,name4)
            subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size13/Combined_Growth.R', location_of_call + '/' , models[0], models[1],models[2],models[3],models[4],models[5],models[6],models[7],models[8],models[9],models[10],models[11],models[12], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11,name12,name13])
            # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4

            ln_num = 0
            growth_info = {}

            Community_growth[models[0].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[1].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[2].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[3].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[4].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[5].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[6].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[7].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[8].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[9].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[10].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[11].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[12].split('/')[-1:][0].split('.RDS')[0]] = []

            for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                ln_num +=1
                if ln_num >1:
                    print (",".join(line.split(',')[:-1][1:]))
                    species = line.split(',')[1].replace('"','')
                    growth = int(line.replace('\n','').split(',')[3])
                    Community_growth[species].append(growth) 

            #subprocess.run(['rm ' + location_of_call + '/' + 'statistics.csv'])
            with open(community_outputfile, 'w') as outputting:
                outputting.write(community_title)
                species_1 = models[0].split('/')[-1:][0].split('.RDS')[0]
                species_2 = models[1].split('/')[-1:][0].split('.RDS')[0]
                species_3 = models[2].split('/')[-1:][0].split('.RDS')[0]
                species_4 = models[3].split('/')[-1:][0].split('.RDS')[0]
                species_5 = models[4].split('/')[-1:][0].split('.RDS')[0]
                species_6 = models[5].split('/')[-1:][0].split('.RDS')[0]
                species_7 = models[6].split('/')[-1:][0].split('.RDS')[0]
                species_8 = models[7].split('/')[-1:][0].split('.RDS')[0]
                species_9 = models[8].split('/')[-1:][0].split('.RDS')[0]
                species_10 = models[9].split('/')[-1:][0].split('.RDS')[0]
                species_11 = models[10].split('/')[-1:][0].split('.RDS')[0]
                species_12 = models[11].split('/')[-1:][0].split('.RDS')[0]
                species_13 = models[12].split('/')[-1:][0].split('.RDS')[0]


                for species in [species_1,species_2,species_3,species_4,species_5,species_6,species_7,species_8,species_9,species_10,species_11,species_12,species_13]:
                    single1 = individual_growths[species][-1:][0]
                    paired1 = Community_growth[species][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    outputting.write(species + '\t' + str(perc1) + '%\n')


            ## Store the data

            all_data_stored[tuple(models)] = [individual_growths,Paired_growths,Community_growth]

            #print ('Completed consortia; ', str(len(all_data_stored.keys())))


        # In[109]:


        # Fourteen species


        # In[ ]:


        if consortia_size_wanted == 14:
            for listi,individual in enumerate(models):
                tmp_models = list(models)
                tmp_models.pop(listi)

                name1 = models[listi].split('/')[-1:][0].split('.RDS')[0]
                name2 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[6].split('/')[-1:][0].split('.RDS')[0]
                name9 = tmp_models[7].split('/')[-1:][0].split('.RDS')[0]
                name10 = tmp_models[8].split('/')[-1:][0].split('.RDS')[0]
                name11 = tmp_models[9].split('/')[-1:][0].split('.RDS')[0]
                name12 = tmp_models[10].split('/')[-1:][0].split('.RDS')[0]
                name13 = tmp_models[11].split('/')[-1:][0].split('.RDS')[0]
                name14 = tmp_models[12].split('/')[-1:][0].split('.RDS')[0]

                print ('Currently growing; ', name1)


                # !!!!!!!!!!!!! Here, the species being studied is given first, then the other ones, then all the names
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size14/Single_Growth.R', location_of_call + '/' , models[listi], tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5],tmp_models[6],tmp_models[7],tmp_models[8],tmp_models[9],tmp_models[10],tmp_models[11],tmp_models[12], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11,name12,name13,name14])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                growth_curve = []
                ln_num = 0
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1]
                        growth_curve.append(int(line.replace('\n','').split(',')[3]))
                individual_growths[name1] = growth_curve

            subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])

            print ('Individual growth curves generated for ;' + str(len(individual_growths)))

                
            all_consortia_combinations = list(itertools.combinations(models, 2))

            Paired_growths = {}

            for comb in all_consortia_combinations:
                tmp_models = list(models[:])
                tmp_models.remove(list(comb)[0])
                tmp_models.remove(list(comb)[1])

                name1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                name2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                name9 = tmp_models[6].split('/')[-1:][0].split('.RDS')[0]
                name10 = tmp_models[7].split('/')[-1:][0].split('.RDS')[0]
                name11 = tmp_models[8].split('/')[-1:][0].split('.RDS')[0]
                name12 = tmp_models[9].split('/')[-1:][0].split('.RDS')[0]
                name13 = tmp_models[10].split('/')[-1:][0].split('.RDS')[0]
                name14 = tmp_models[11].split('/')[-1:][0].split('.RDS')[0]
                print ('Strains being grown together;')
                print( list(comb)[0].replace('"','').split('/')[-1:][0].split('.RDS')[0], list(comb)[1].replace('"','').split('/')[-1:][0].split('.RDS')[0], '; names=(', name1,',',name2,')') # Only highlighting those species being paired
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size14/Paired_Growth.R', location_of_call + '/' , list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5],tmp_models[6],tmp_models[7],tmp_models[8],tmp_models[9],tmp_models[10],tmp_models[11], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11,name12,name13,name14])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                ln_num = 0
                growth_info = {}
                growth_info[list(comb)[0].split('/')[-1:][0].split('.RDS')[0]] = []
                growth_info[list(comb)[1].split('/')[-1:][0].split('.RDS')[0]] = []
                for line in open(location_of_call + '/' + 'statistics.csv'):
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1].replace('"','').split('/')[-1:][0].split('.RDS')[0]
                        growth = int(line.replace('\n','').split(',')[3])
                        growth_info[species].append(growth) 
                Paired_growths[comb] = growth_info

            subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])

            ### Identify paired interactions
            with open(paired_outputfile, 'w') as outputting:
                outputting.write(paired_title)
                for comb, data in Paired_growths.items():
                    print ('Calculating interaction between; ')
                    #print (comb)
                    species_1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                    species_2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                    print (species_1)
                    print (species_2)

                    # Calculate interaction type for species1

                    single1 = individual_growths[species_1][-1:][0]
                    paired1 = data[species_1][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    if perc1 >10.0:
                        interaction1 = 'positive'
                    elif perc1 < -10.0:
                        interaction1 = 'negative'
                    else:
                        interaction1 = 'none'


                    # Calculate interaction type for species2

                    single2 = individual_growths[species_2][-1:][0]
                    paired2 = data[species_2][-1:][0]

                    perc2 = ((single2/paired2)*100)-100
                    print (perc1, perc2)

                    if perc2 >10.0:
                        interaction2 = 'positive'
                    elif perc2 < -10.0:
                        interaction2 = 'negative'
                    else:
                        interaction2 = 'none'


                    # Overall interaction type

                    if interaction1 == 'positive':
                        if interaction2 == 'positive':
                            overall_interaction = 'mutualism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'none':
                            overall_interaction = 'commensalism'

                    elif interaction1 == 'negative':
                        if interaction2 == 'positive':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'competition'
                        elif interaction2 == 'none':
                            overall_interaction = 'amensalism'


                    elif interaction1 == 'none':
                        if interaction2 == 'positive':
                            overall_interaction = 'commensalism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'amensalism'
                        elif interaction2 == 'none':
                            overall_interaction = 'neutralism'

                    outputting.write(species_1 + ' and ' + species_2 + '\t' + overall_interaction + '\n')
                    print (species_1 + ' and ' + species_2 + '; ' + overall_interaction + '\n')
            
            print ('Running whole community model.')
            Community_growth = {}

            name1 = models[0].split('/')[-1:][0].split('.RDS')[0]
            name2 = models[1].split('/')[-1:][0].split('.RDS')[0]
            name3 = models[2].split('/')[-1:][0].split('.RDS')[0]
            name4 = models[3].split('/')[-1:][0].split('.RDS')[0]
            name5 = models[4].split('/')[-1:][0].split('.RDS')[0]
            name6 = models[5].split('/')[-1:][0].split('.RDS')[0]
            name7 = models[6].split('/')[-1:][0].split('.RDS')[0]
            name8 = models[7].split('/')[-1:][0].split('.RDS')[0]
            name9 = models[8].split('/')[-1:][0].split('.RDS')[0]
            name10 = models[9].split('/')[-1:][0].split('.RDS')[0]
            name11 = models[10].split('/')[-1:][0].split('.RDS')[0]
            name12 = models[11].split('/')[-1:][0].split('.RDS')[0]
            name13 = models[12].split('/')[-1:][0].split('.RDS')[0]
            name14 = models[13].split('/')[-1:][0].split('.RDS')[0]
            try:
                subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])
            except FileNotFoundError:
                lol = 0

            #print( list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1], name1,name2,name3,name4)
            subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size14/Combined_Growth.R', location_of_call + '/' , models[0], models[1],models[2],models[3],models[4],models[5],models[6],models[7],models[8],models[9],models[10],models[11],models[12],models[13], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11,name12,name13,name14])
            # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4

            ln_num = 0
            growth_info = {}

            Community_growth[models[0].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[1].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[2].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[3].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[4].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[5].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[6].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[7].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[8].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[9].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[10].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[11].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[12].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[13].split('/')[-1:][0].split('.RDS')[0]] = []

            for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                ln_num +=1
                if ln_num >1:
                    print (",".join(line.split(',')[:-1][1:]))
                    species = line.split(',')[1].replace('"','')
                    growth = int(line.replace('\n','').split(',')[3])
                    Community_growth[species].append(growth) 

            #subprocess.run(['rm ' + location_of_call + '/' + 'statistics.csv'])
            with open(community_outputfile, 'w') as outputting:
                outputting.write(community_title)

                species_1 = models[0].split('/')[-1:][0].split('.RDS')[0]
                species_2 = models[1].split('/')[-1:][0].split('.RDS')[0]
                species_3 = models[2].split('/')[-1:][0].split('.RDS')[0]
                species_4 = models[3].split('/')[-1:][0].split('.RDS')[0]
                species_5 = models[4].split('/')[-1:][0].split('.RDS')[0]
                species_6 = models[5].split('/')[-1:][0].split('.RDS')[0]
                species_7 = models[6].split('/')[-1:][0].split('.RDS')[0]
                species_8 = models[7].split('/')[-1:][0].split('.RDS')[0]
                species_9 = models[8].split('/')[-1:][0].split('.RDS')[0]
                species_10 = models[9].split('/')[-1:][0].split('.RDS')[0]
                species_11 = models[10].split('/')[-1:][0].split('.RDS')[0]
                species_12 = models[11].split('/')[-1:][0].split('.RDS')[0]
                species_13 = models[12].split('/')[-1:][0].split('.RDS')[0]
                species_14 = models[13].split('/')[-1:][0].split('.RDS')[0]


                for species in [species_1,species_2,species_3,species_4,species_5,species_6,species_7,species_8,species_9,species_10,species_11,species_12,species_13,species_14]:
                    single1 = individual_growths[species][-1:][0]
                    paired1 = Community_growth[species][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    outputting.write(species + '\t' + str(perc1) + '%\n')


            ## Store the data

            all_data_stored[tuple(models)] = [individual_growths,Paired_growths,Community_growth]

            #print ('Completed consortia; ', str(len(all_data_stored.keys())))


        # In[110]:


        # Fifthteen species


        # In[17]:


        if consortia_size_wanted == 15:
            for listi,individual in enumerate(models):
                tmp_models = list(models)
                tmp_models.pop(listi)

                name1 = models[listi].split('/')[-1:][0].split('.RDS')[0]
                name2 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[6].split('/')[-1:][0].split('.RDS')[0]
                name9 = tmp_models[7].split('/')[-1:][0].split('.RDS')[0]
                name10 = tmp_models[8].split('/')[-1:][0].split('.RDS')[0]
                name11 = tmp_models[9].split('/')[-1:][0].split('.RDS')[0]
                name12 = tmp_models[10].split('/')[-1:][0].split('.RDS')[0]
                name13 = tmp_models[11].split('/')[-1:][0].split('.RDS')[0]
                name14 = tmp_models[12].split('/')[-1:][0].split('.RDS')[0]
                name15 = tmp_models[13].split('/')[-1:][0].split('.RDS')[0]

                print ('Currently growing; ', name1)


                # !!!!!!!!!!!!! Here, the species being studied is given first, then the other ones, then all the names
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size15/Single_Growth.R', location_of_call + '/' , models[listi], tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5],tmp_models[6],tmp_models[7],tmp_models[8],tmp_models[9],tmp_models[10],tmp_models[11],tmp_models[12],tmp_models[13], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11,name12,name13,name14,name15])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                growth_curve = []
                ln_num = 0
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1]
                        growth_curve.append(int(line.replace('\n','').split(',')[3]))
                individual_growths[name1] = growth_curve

            subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])

            print ('Individual growth curves generated for ;' + str(len(individual_growths)))

                
            all_consortia_combinations = list(itertools.combinations(models, 2))

            Paired_growths = {}

            for comb in all_consortia_combinations:
                tmp_models = list(models[:])
                tmp_models.remove(list(comb)[0])
                tmp_models.remove(list(comb)[1])

                name1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                name2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                name3 = tmp_models[0].split('/')[-1:][0].split('.RDS')[0]
                name4 = tmp_models[1].split('/')[-1:][0].split('.RDS')[0]
                name5 = tmp_models[2].split('/')[-1:][0].split('.RDS')[0]
                name6 = tmp_models[3].split('/')[-1:][0].split('.RDS')[0]
                name7 = tmp_models[4].split('/')[-1:][0].split('.RDS')[0]
                name8 = tmp_models[5].split('/')[-1:][0].split('.RDS')[0]
                name9 = tmp_models[6].split('/')[-1:][0].split('.RDS')[0]
                name10 = tmp_models[7].split('/')[-1:][0].split('.RDS')[0]
                name11 = tmp_models[8].split('/')[-1:][0].split('.RDS')[0]
                name12 = tmp_models[9].split('/')[-1:][0].split('.RDS')[0]
                name13 = tmp_models[10].split('/')[-1:][0].split('.RDS')[0]
                name14 = tmp_models[11].split('/')[-1:][0].split('.RDS')[0]
                name15 = tmp_models[12].split('/')[-1:][0].split('.RDS')[0]
                print ('Strains being grown together;')
                print( list(comb)[0].replace('"','').split('/')[-1:][0].split('.RDS')[0], list(comb)[1].replace('"','').split('/')[-1:][0].split('.RDS')[0], '; names=(', name1,',',name2,')') # Only highlighting those species being paired
                subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size15/Paired_Growth.R', location_of_call + '/' , list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1],tmp_models[2],tmp_models[3],tmp_models[4],tmp_models[5],tmp_models[6],tmp_models[7],tmp_models[8],tmp_models[9],tmp_models[10],tmp_models[11],tmp_models[12], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11,name12,name13,name14,name15])
                # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4
                ln_num = 0
                growth_info = {}
                growth_info[list(comb)[0].split('/')[-1:][0].split('.RDS')[0]] = []
                growth_info[list(comb)[1].split('/')[-1:][0].split('.RDS')[0]] = []
                for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                    ln_num +=1
                    if ln_num >1:
                        print (",".join(line.split(',')[:-1][1:]))
                        species = line.split(',')[1].replace('"','').split('/')[-1:][0].split('.RDS')[0]
                        growth = int(line.replace('\n','').split(',')[3])
                        growth_info[species].append(growth) 
                Paired_growths[comb] = growth_info

            subprocess.run(['rm', location_of_call + '/' + 'statistics.csv'])

            ### Identify paired interactions
            with open(paired_outputfile, 'w') as outputting:
                for comb, data in Paired_growths.items():
                    print ('Calculating interaction between; ')
                    #print (comb)
                    species_1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
                    species_2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]
                    print (species_1)
                    print (species_2)

                    # Calculate interaction type for species1

                    single1 = individual_growths[species_1][-1:][0]
                    paired1 = data[species_1][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    if perc1 >10.0:
                        interaction1 = 'positive'
                    elif perc1 < -10.0:
                        interaction1 = 'negative'
                    else:
                        interaction1 = 'none'


                    # Calculate interaction type for species2

                    single2 = individual_growths[species_2][-1:][0]
                    paired2 = data[species_2][-1:][0]

                    perc2 = ((single2/paired2)*100)-100
                    print (perc1, perc2)

                    if perc2 >10.0:
                        interaction2 = 'positive'
                    elif perc2 < -10.0:
                        interaction2 = 'negative'
                    else:
                        interaction2 = 'none'


                    # Overall interaction type

                    if interaction1 == 'positive':
                        if interaction2 == 'positive':
                            overall_interaction = 'mutualism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'none':
                            overall_interaction = 'commensalism'

                    elif interaction1 == 'negative':
                        if interaction2 == 'positive':
                            overall_interaction = 'parasitism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'competition'
                        elif interaction2 == 'none':
                            overall_interaction = 'amensalism'


                    elif interaction1 == 'none':
                        if interaction2 == 'positive':
                            overall_interaction = 'commensalism'
                        elif interaction2 == 'negative':
                            overall_interaction = 'amensalism'
                        elif interaction2 == 'none':
                            overall_interaction = 'neutralism'

                    outputting.write(species_1 + ' and ' + species_2 + '; ' + overall_interaction + '\n')
                    print (species_1 + ' and ' + species_2 + '; ' + overall_interaction + '\n')
            
            print ('Running whole community model.')
            Community_growth = {}

            name1 = models[0].split('/')[-1:][0].split('.RDS')[0]
            name2 = models[1].split('/')[-1:][0].split('.RDS')[0]
            name3 = models[2].split('/')[-1:][0].split('.RDS')[0]
            name4 = models[3].split('/')[-1:][0].split('.RDS')[0]
            name5 = models[4].split('/')[-1:][0].split('.RDS')[0]
            name6 = models[5].split('/')[-1:][0].split('.RDS')[0]
            name7 = models[6].split('/')[-1:][0].split('.RDS')[0]
            name8 = models[7].split('/')[-1:][0].split('.RDS')[0]
            name9 = models[8].split('/')[-1:][0].split('.RDS')[0]
            name10 = models[9].split('/')[-1:][0].split('.RDS')[0]
            name11 = models[10].split('/')[-1:][0].split('.RDS')[0]
            name12 = models[11].split('/')[-1:][0].split('.RDS')[0]
            name13 = models[12].split('/')[-1:][0].split('.RDS')[0]
            name14 = models[13].split('/')[-1:][0].split('.RDS')[0]
            name15 = models[14].split('/')[-1:][0].split('.RDS')[0]
            try:
                subprocess.run(['rm',  location_of_call + '/' + 'statistics.csv'])
            except FileNotFoundError:
                lol = 0

            #print( list(comb)[0], list(comb)[1],tmp_models[0],tmp_models[1], name1,name2,name3,name4)
            subprocess.run(['Rscript', location_of_file + '/' + 'scripts/Modelling/Size15/Combined_Growth.R', location_of_call + '/' , models[0], models[1],models[2],models[3],models[4],models[5],models[6],models[7],models[8],models[9],models[10],models[11],models[12],models[13],models[14], name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11,name12,name13,name14,name15])
            # Commands must be Rscript Single_Growth.R working_directory model1 model2 model3 model4 Name1 Name2 Name3 Name4

            ln_num = 0
            growth_info = {}

            Community_growth[models[0].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[1].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[2].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[3].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[4].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[5].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[6].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[7].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[8].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[9].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[10].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[11].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[12].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[13].split('/')[-1:][0].split('.RDS')[0]] = []
            Community_growth[models[14].split('/')[-1:][0].split('.RDS')[0]] = []

            for line in open(location_of_call + '/' + 'statistics.csv','r'): #change
                ln_num +=1
                if ln_num >1:
                    print (",".join(line.split(',')[:-1][1:]))
                    species = line.split(',')[1].replace('"','')
                    growth = int(line.replace('\n','').split(',')[3])
                    Community_growth[species].append(growth) 

            #subprocess.run(['rm ' + location_of_call + '/' + 'statistics.csv'])
            with open(community_outputfile, 'w') as outputting:
                outputting.write(community_title)

                species_1 = models[0].split('/')[-1:][0].split('.RDS')[0]
                species_2 = models[1].split('/')[-1:][0].split('.RDS')[0]
                species_3 = models[2].split('/')[-1:][0].split('.RDS')[0]
                species_4 = models[3].split('/')[-1:][0].split('.RDS')[0]
                species_5 = models[4].split('/')[-1:][0].split('.RDS')[0]
                species_6 = models[5].split('/')[-1:][0].split('.RDS')[0]
                species_7 = models[6].split('/')[-1:][0].split('.RDS')[0]
                species_8 = models[7].split('/')[-1:][0].split('.RDS')[0]
                species_9 = models[8].split('/')[-1:][0].split('.RDS')[0]
                species_10 = models[9].split('/')[-1:][0].split('.RDS')[0]
                species_11 = models[10].split('/')[-1:][0].split('.RDS')[0]
                species_12 = models[11].split('/')[-1:][0].split('.RDS')[0]
                species_13 = models[12].split('/')[-1:][0].split('.RDS')[0]
                species_14 = models[13].split('/')[-1:][0].split('.RDS')[0]
                species_15 = models[14].split('/')[-1:][0].split('.RDS')[0]


                for species in [species_1,species_2,species_3,species_4,species_5,species_6,species_7,species_8,species_9,species_10,species_11,species_12,species_13,species_14,species_15]:
                    single1 = individual_growths[species][-1:][0]
                    paired1 = Community_growth[species][-1:][0]

                    perc1 = ((single1/paired1)*100)-100

                    outputting.write(species + '\t' + str(perc1) + '%\n')


            ## Store the data

            all_data_stored[tuple(models)] = [individual_growths,Paired_growths,Community_growth]

            #print ('Completed consortia; ', str(len(all_data_stored.keys())))


        # # Check is models were sufficient

        # In[113]:


        print ('Assessing paired interactions....')
        total_interactions = 0
        not_negative_interactions = 0

        species_not_negative_interactions = {}
        for species in individual_growths.keys():
            species_not_negative_interactions[species] = 0

        for comb, data in Paired_growths.items():
            #print (comb)
            species_1 = list(comb)[0].split('/')[-1:][0].split('.RDS')[0]
            species_2 = list(comb)[1].split('/')[-1:][0].split('.RDS')[0]

            # Calculate interaction type for species1

            single1 = individual_growths[species_1][-1:][0]
            paired1 = data[species_1][-1:][0]

            perc1 = ((single1/paired1)*100)-100
            
            total_interactions +=2
            
            if perc1 >10.0:
                interaction1 = 'positive'
                not_negative_interactions +=1
                species_not_negative_interactions[species_1] +=1
            elif perc1 < -10.0:
                interaction1 = 'negative'
            else:
                interaction1 = 'none'
                not_negative_interactions +=1
                species_not_negative_interactions[species_1] +=1


            # Calculate interaction type for species2

            single2 = individual_growths[species_2][-1:][0]
            paired2 = data[species_2][-1:][0]

            perc2 = ((single2/paired2)*100)-100

            if perc2 >10.0:
                interaction2 = 'positive'
                not_negative_interactions +=1
                species_not_negative_interactions[species_2] +=1
            elif perc2 < -10.0:
                interaction2 = 'negative'
            else:
                interaction2 = 'none'
                not_negative_interactions +=1
                species_not_negative_interactions[species_2] +=1


        percentage_of_interactions_not_negative = (not_negative_interactions/total_interactions)*100


        ### Determine if the community is OK based on paired growth

        if percentage_of_interactions_not_negative > 50.0:
            print ('Community has passed phase 1 stability check; the majority of paired interactions are not negative.')
            if all(value != 0 for value in species_not_negative_interactions.values()): # All species have one non-negative interaction
                print ('Community has passed phase 2 stability check; all isolates have atleast one non-negative interaction.')
                passed_phase_one_stability = True
            else:
                print ('Community is unstable; some isolates have only negative interactions!')
                passed_phase_one_stability = False
        else:
            print ('Community is unstable: over 50% of paired interactions are negative!')
            passed_phase_one_stability = False

        if passed_phase_one_stability == True:    
            did_all_grow_in_community = 0
            for species, growth_data in Community_growth.items():
                if growth_data[-1:][0] > 10: # Is the value after growth greater than the starting number of cells?
                        did_all_grow_in_community +=1

            if did_all_grow_in_community == len(models):
                print ('Community has passed phase 3 stability check; All isolates were observed to grow in the community model.')
                passed_phase_two_stability = True
            else:
                print ('Community is unstable; not all isolates were observed to grow in the community model.')
                passed_phase_two_stability = False

        if passed_phase_one_stability == False or passed_phase_two_stability == False: # The community predicted isnt good, so another round is needed
            modelling_iteration += 1
            subprocess.run(['mkdir', location_of_call + '/'+ output_prefix + 'Modelling-iteration-' + str(modelling_iteration)]) ###!!!CHANGE!!!
            subprocess.run(['mv',location_of_call + '/'+ output_prefix +'-Paired_Interactions.tsv',location_of_call + '/'+ output_prefix + 'Modelling-iteration-' + str(modelling_iteration)])
            subprocess.run(['mv',location_of_call + '/'+ output_prefix +'-Community-wide_response.tsv',location_of_call + '/'+ output_prefix + 'Modelling-iteration-' + str(modelling_iteration)])
            consortia_mean_Scores.pop(max_match)


print ('::::::: Step 4 complete: Metabolic modelling and interaction analysis.')

print (':::::::: Outputting analysis to user.')

# Final step to generate figures and summary files

tax_outputfile = output_prefix +'-SynCom_Taxonomy.tsv'
with open(tax_outputfile, 'w') as outputting_taxonomy:
    outputting_taxonomy.write('Strain\tTaxonomic assignment\n')
    if taxonomic_filtering == True:
        print ('The taxonomic assignments of the selected consortia members are; ')
        for i in consortia_selected:
            print (i +' = ' + taxonomy.loc[i]['classification'])
            outputting_taxonomy.write(i + '\t' + taxonomy.loc[i]['classification'] + '\n')


def sample_consortia_accounted(column, consortia):  
    sample_data = samples[column]
    samples_pfams = list(sample_data.loc[sample_data == 1].index)
    consortia_pfams = []
    for genome in consortia:   
        genome_data = genomes[genome] # This section creates the single genome file
        redundant_functional_list = consortia_pfams + list(genome_data.loc[genome_data == 1].index)
        consortia_pfams = list(set(redundant_functional_list))


    match = 0
    mismatch = 0

    match = len(list(set(samples_pfams).intersection(consortia_pfams))) # Identifies the common 
    mismatch = len(consortia_pfams)-match
    account_perc = (match / len(samples_pfams)) * 100 
    return match, mismatch, account_perc


pairwise_sample_accounted = {}
pairwise_sample_mismatches = {}

for column, v in Stored_first_round_consortia.items():
    consortia_matches, consortia_mismatches, cul_accounted,cul_mismatches, consortia = v
    
    cut_off_consortia = consortia[0:consortia_size_wanted]
    #print (cut_off_consortia)
    #print (column)
    reduced_match, reduced_mimatch, reduced_accum = sample_consortia_accounted(column, cut_off_consortia)

    final_match, final_mimatch, final_accum = sample_consortia_accounted(column, consortia_selected) 
    

    pairwise_sample_accounted[column] = [cul_accounted[-1:][0], reduced_accum, final_accum]

    
    pairwise_sample_mismatches[column] = [cul_mismatches[-1:][0], reduced_mimatch, final_mimatch]

    
plotting_data = pd.DataFrame.from_dict(pairwise_sample_accounted, orient='index', columns=['Sample specific (size = 20)', 'Sample specific (size =' + str(consortia_size_wanted) + ')', 'Group-wide selection'])


plt.figure(figsize=(8, 5))
# Create violin plots without mini-boxplots inside.
ax = sns.violinplot( data=plotting_data,
                    palette=['#d7eaf3', '#77b5d9', '#14397d'], 
                    cut=0, inner=None, linewidth=0.03,color='k')
# Clip the right half of each violin.
for item in ax.collections:
    x0, y0, width, height = item.get_paths()[0].get_extents().bounds
    item.set_clip_path(plt.Rectangle((x0, y0), width/2, height,
                       transform=ax.transData))
# Create strip plots with partially transparent points of different colors depending on the group.
num_items = len(ax.collections)
sns.stripplot(  data=plotting_data,
              palette=['#d7eaf3', '#77b5d9', '#14397d'], alpha=0.4, size=7)
# Shift each strip plot strictly below the correponding volin.
for item in ax.collections[num_items:]:
    item.set_offsets(item.get_offsets() + 0.15)
# Create narrow boxplots on top of the corresponding violin and strip plots, with thick lines, the mean values, without the outliers.
sns.boxplot( data=plotting_data, width=0.25,
            showfliers=False, showmeans=True, 
            meanprops=dict(marker='o', markerfacecolor='darkred',
                           markersize=10, zorder=3),
            boxprops=dict(facecolor=(0,0,0,0), 
                          linewidth=1, zorder=3),
            whiskerprops=dict(linewidth=1),
            capprops=dict(linewidth=1),
            medianprops=dict(linewidth=1))



ax.set( ylabel='Accounted Pfams (%)')
#plt.show()
plt.savefig(location_of_call + '/'+ output_prefix + '-Accounted_Pfams_across_consortia.pdf', bbox_inches='tight')
plt.close()
plt.figure().clear()

plotting_data = pd.DataFrame.from_dict(pairwise_sample_mismatches, orient='index', columns=['Sample specific (size = 20)', 'Sample specific (size =' + str(consortia_size_wanted) + ')', 'Group-wide selection'])

plt.figure(figsize=(8, 5))
# Create violin plots without mini-boxplots inside.
ax = sns.violinplot( data=plotting_data,
                    palette=['#d7eaf3', '#77b5d9', '#14397d'], 
                    cut=0, inner=None, linewidth=0.03,color='k')
# Clip the right half of each violin.
for item in ax.collections:
    x0, y0, width, height = item.get_paths()[0].get_extents().bounds
    item.set_clip_path(plt.Rectangle((x0, y0), width/2, height,
                       transform=ax.transData))
# Create strip plots with partially transparent points of different colors depending on the group.
num_items = len(ax.collections)
sns.stripplot(  data=plotting_data,
              palette=['#d7eaf3', '#77b5d9', '#14397d'], alpha=0.4, size=7)
# Shift each strip plot strictly below the correponding volin.
for item in ax.collections[num_items:]:
    item.set_offsets(item.get_offsets() + 0.15)
# Create narrow boxplots on top of the corresponding violin and strip plots, with thick lines, the mean values, without the outliers.
sns.boxplot( data=plotting_data, width=0.25,
            showfliers=False, showmeans=True, 
            meanprops=dict(marker='o', markerfacecolor='darkred',
                           markersize=10, zorder=3),
            boxprops=dict(facecolor=(0,0,0,0), 
                          linewidth=1, zorder=3),
            whiskerprops=dict(linewidth=1),
            capprops=dict(linewidth=1),
            medianprops=dict(linewidth=1))



ax.set( ylabel='Mismatching Pfams')
#plt.show()

plt.savefig(location_of_call + '/'+ output_prefix + '-Mismatches_across_consortia.pdf', bbox_inches='tight')
plt.close()
plt.figure().clear()

#  requires statistics.csv generated during Step 4

if models_folder is not None and models_folder != " ":
    df = pd.read_csv('statistics.csv')
    species = list(set(list(df['data.species'])))
    times = list(set(list(df['data.time'])))
    raw_data = {}

    for sp in species:
        growth = []
        for gtime in times:
            growth.append(int(df.loc[df['data.species'] == sp].loc[df['data.time'] == gtime]['data.value']))
        raw_data[sp] = growth

    r = times
    df2 = pd.DataFrame(raw_data)
    
        
        
    if consortia_size_wanted == 2:
        # From raw value to percentage
        totals = [i+j for i,j in zip(df2[species[0]], df2[species[1]])]
        greenBars = [i / j * 100 for i,j in zip(df2[species[0]], totals)]
        orangeBars = [i / j * 100 for i,j in zip(df2[species[1]], totals)]

        # plot
        barWidth = 0.85
        # Create green Bars
        plt.bar(r, greenBars, edgecolor='white', width=barWidth, label = species[0])
        # Create orange Bars
        plt.bar(r, orangeBars, bottom=greenBars, edgecolor='white', width=barWidth, label = species[1])
        # Create blue Bars

        # Custom x axis
        plt.xlabel("Time (h)")

        plt.ylabel('Relative abundance (%)')
        plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)

        # Show graphic
        plt.savefig(location_of_call + '/'+ output_prefix + '-Relative_Abundance_growth.pdf', bbox_inches='tight')
        plt.close()
        plt.figure().clear()
        
        
    if consortia_size_wanted == 3:
        # From raw value to percentage
        totals = [i+j+k for i,j,k in zip(df2[species[0]], df2[species[1]], df2[species[2]])]
        greenBars = [i / j * 100 for i,j in zip(df2[species[0]], totals)]
        orangeBars = [i / j * 100 for i,j in zip(df2[species[1]], totals)]
        threeBars = [i / j * 100 for i,j in zip(df2[species[2]], totals)]

        # plot
        barWidth = 0.85
        # Create green Bars
        plt.bar(r, greenBars, edgecolor='white', width=barWidth, label = species[0])
        # Create orange Bars
        plt.bar(r, orangeBars, bottom=greenBars, edgecolor='white', width=barWidth, label = species[1])
        # Create three Bars
        plt.bar(r, threeBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], edgecolor='white', width=barWidth, label = species[2])
        

        # Custom x axis
        plt.xlabel("Time (h)")

        plt.ylabel('Relative abundance (%)')
        plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)

        # Show graphic
        plt.savefig(location_of_call + '/'+ output_prefix + '-Relative_Abundance_growth.pdf', bbox_inches='tight')
        plt.close()
        plt.figure().clear()

    if consortia_size_wanted == 4:
        # From raw value to percentage
        totals = [i+j+k+l for i,j,k,l in zip(df2[species[0]], df2[species[1]], df2[species[2]], df2[species[3]])]
        greenBars = [i / j * 100 for i,j in zip(df2[species[0]], totals)]
        orangeBars = [i / j * 100 for i,j in zip(df2[species[1]], totals)]
        threeBars = [i / j * 100 for i,j in zip(df2[species[2]], totals)]
        fourBars = [i / j * 100 for i,j in zip(df2[species[3]], totals)]

        # plot
        barWidth = 0.85
        # Create green Bars
        plt.bar(r, greenBars, edgecolor='white', width=barWidth, label = species[0])
        # Create orange Bars
        plt.bar(r, orangeBars, bottom=greenBars, edgecolor='white', width=barWidth, label = species[1])
        # Create three Bars
        plt.bar(r, threeBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], edgecolor='white', width=barWidth, label = species[2])
        # Create four Bars
        plt.bar(r, fourBars, bottom=[i+j+k for i,j,k in zip(greenBars, orangeBars,threeBars)], edgecolor='white', width=barWidth, label = species[3])
        

        # Custom x axis
        plt.xlabel("Time (h)")

        plt.ylabel('Relative abundance (%)')
        plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)

        # Show graphic
        plt.savefig(location_of_call + '/'+ output_prefix + '-Relative_Abundance_growth.pdf', bbox_inches='tight')
        plt.close()
        plt.figure().clear()
        
        
    if consortia_size_wanted == 5:
        # From raw value to percentage
        totals = [i+j+k+l+m for i,j,k,l,m in zip(df2[species[0]], df2[species[1]], df2[species[2]], df2[species[3]], df2[species[4]])]
        greenBars = [i / j * 100 for i,j in zip(df2[species[0]], totals)]
        orangeBars = [i / j * 100 for i,j in zip(df2[species[1]], totals)]
        threeBars = [i / j * 100 for i,j in zip(df2[species[2]], totals)]
        fourBars = [i / j * 100 for i,j in zip(df2[species[3]], totals)]
        fiveBars = [i / j * 100 for i,j in zip(df2[species[4]], totals)]

        # plot
        barWidth = 0.85
        # Create green Bars
        plt.bar(r, greenBars, edgecolor='white', width=barWidth, label = species[0])
        # Create orange Bars
        plt.bar(r, orangeBars, bottom=greenBars, edgecolor='white', width=barWidth, label = species[1])
        # Create three Bars
        plt.bar(r, threeBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], edgecolor='white', width=barWidth, label = species[2])
        # Create four Bars
        plt.bar(r, fourBars, bottom=[i+j+k for i,j,k in zip(greenBars, orangeBars,threeBars)], edgecolor='white', width=barWidth, label = species[3])
        # Create five Bars
        plt.bar(r, fiveBars, bottom=[i+j+k+l for i,j,k,l in zip(greenBars, orangeBars,threeBars,fourBars)], edgecolor='white', width=barWidth, label = species[4])
        

        # Custom x axis
        plt.xlabel("Time (h)")

        plt.ylabel('Relative abundance (%)')
        plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)

        # Show graphic
        plt.savefig(location_of_call + '/'+ output_prefix + '-Relative_Abundance_growth.pdf', bbox_inches='tight')
        plt.close()
        plt.figure().clear()
        
        
    if consortia_size_wanted == 6:
        # From raw value to percentage
        totals = [i+j+k+l+m+n for i,j,k,l,m,n in zip(df2[species[0]], df2[species[1]], df2[species[2]], df2[species[3]], df2[species[4]], df2[species[5]])]
        greenBars = [i / j * 100 for i,j in zip(df2[species[0]], totals)]
        orangeBars = [i / j * 100 for i,j in zip(df2[species[1]], totals)]
        threeBars = [i / j * 100 for i,j in zip(df2[species[2]], totals)]
        fourBars = [i / j * 100 for i,j in zip(df2[species[3]], totals)]
        fiveBars = [i / j * 100 for i,j in zip(df2[species[4]], totals)]
        sixBars = [i / j * 100 for i,j in zip(df2[species[5]], totals)]

        # plot
        barWidth = 0.85
        # Create green Bars
        plt.bar(r, greenBars, edgecolor='white', width=barWidth, label = species[0])
        # Create orange Bars
        plt.bar(r, orangeBars, bottom=greenBars, edgecolor='white', width=barWidth, label = species[1])
        # Create three Bars
        plt.bar(r, threeBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], edgecolor='white', width=barWidth, label = species[2])
        # Create four Bars
        plt.bar(r, fourBars, bottom=[i+j+k for i,j,k in zip(greenBars, orangeBars,threeBars)], edgecolor='white', width=barWidth, label = species[3])
        # Create five Bars
        plt.bar(r, fiveBars, bottom=[i+j+k+l for i,j,k,l in zip(greenBars, orangeBars,threeBars,fourBars)], edgecolor='white', width=barWidth, label = species[4])
        # Create six Bars
        plt.bar(r, sixBars, bottom=[i+j+k+l+m for i,j,k,l,m in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars)], edgecolor='white', width=barWidth, label = species[5])
        

        # Custom x axis
        plt.xlabel("Time (h)")

        plt.ylabel('Relative abundance (%)')
        plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)

        # Show graphic
        plt.savefig(location_of_call + '/'+ output_prefix + '-Relative_Abundance_growth.pdf', bbox_inches='tight')
        plt.close()
        plt.figure().clear()
        
    if consortia_size_wanted == 7:
        # From raw value to percentage
        totals = [i+j+k+l+m+n+o for i,j,k,l,m,n,o in zip(df2[species[0]], df2[species[1]], df2[species[2]], df2[species[3]], df2[species[4]], df2[species[5]], df2[species[6]])]
        greenBars = [i / j * 100 for i,j in zip(df2[species[0]], totals)]
        orangeBars = [i / j * 100 for i,j in zip(df2[species[1]], totals)]
        threeBars = [i / j * 100 for i,j in zip(df2[species[2]], totals)]
        fourBars = [i / j * 100 for i,j in zip(df2[species[3]], totals)]
        fiveBars = [i / j * 100 for i,j in zip(df2[species[4]], totals)]
        sixBars = [i / j * 100 for i,j in zip(df2[species[5]], totals)]
        sevenBars = [i / j * 100 for i,j in zip(df2[species[6]], totals)]

        # plot
        barWidth = 0.85
        # Create green Bars
        plt.bar(r, greenBars, edgecolor='white', width=barWidth, label = species[0])
        # Create orange Bars
        plt.bar(r, orangeBars, bottom=greenBars, edgecolor='white', width=barWidth, label = species[1])
        # Create three Bars
        plt.bar(r, threeBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], edgecolor='white', width=barWidth, label = species[2])
        # Create four Bars
        plt.bar(r, fourBars, bottom=[i+j+k for i,j,k in zip(greenBars, orangeBars,threeBars)], edgecolor='white', width=barWidth, label = species[3])
        # Create five Bars
        plt.bar(r, fiveBars, bottom=[i+j+k+l for i,j,k,l in zip(greenBars, orangeBars,threeBars,fourBars)], edgecolor='white', width=barWidth, label = species[4])
        # Create six Bars
        plt.bar(r, sixBars, bottom=[i+j+k+l+m for i,j,k,l,m in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars)], edgecolor='white', width=barWidth, label = species[5])
        # Create seven Bars
        plt.bar(r, sevenBars, bottom=[i+j+k+l+m+n for i,j,k,l,m,n in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars)], edgecolor='white', width=barWidth, label = species[6])
        

        # Custom x axis
        plt.xlabel("Time (h)")

        plt.ylabel('Relative abundance (%)')
        plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)

        # Show graphic
        plt.savefig(location_of_call + '/'+ output_prefix + '-Relative_Abundance_growth.pdf', bbox_inches='tight')
        plt.close()
        plt.figure().clear()
        
        
        
    if consortia_size_wanted == 8:
        # From raw value to percentage
        totals = [i+j+k+l+m+n+o+p for i,j,k,l,m,n,o,p in zip(df2[species[0]], df2[species[1]], df2[species[2]], df2[species[3]], df2[species[4]], df2[species[5]], df2[species[6]], df2[species[7]])]
        greenBars = [i / j * 100 for i,j in zip(df2[species[0]], totals)]
        orangeBars = [i / j * 100 for i,j in zip(df2[species[1]], totals)]
        threeBars = [i / j * 100 for i,j in zip(df2[species[2]], totals)]
        fourBars = [i / j * 100 for i,j in zip(df2[species[3]], totals)]
        fiveBars = [i / j * 100 for i,j in zip(df2[species[4]], totals)]
        sixBars = [i / j * 100 for i,j in zip(df2[species[5]], totals)]
        sevenBars = [i / j * 100 for i,j in zip(df2[species[6]], totals)]
        eightBars = [i / j * 100 for i,j in zip(df2[species[7]], totals)]

        # plot
        barWidth = 0.85
        # Create green Bars
        plt.bar(r, greenBars, edgecolor='white', width=barWidth, label = species[0])
        # Create orange Bars
        plt.bar(r, orangeBars, bottom=greenBars, edgecolor='white', width=barWidth, label = species[1])
        # Create three Bars
        plt.bar(r, threeBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], edgecolor='white', width=barWidth, label = species[2])
        # Create four Bars
        plt.bar(r, fourBars, bottom=[i+j+k for i,j,k in zip(greenBars, orangeBars,threeBars)], edgecolor='white', width=barWidth, label = species[3])
        # Create five Bars
        plt.bar(r, fiveBars, bottom=[i+j+k+l for i,j,k,l in zip(greenBars, orangeBars,threeBars,fourBars)], edgecolor='white', width=barWidth, label = species[4])
        # Create six Bars
        plt.bar(r, sixBars, bottom=[i+j+k+l+m for i,j,k,l,m in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars)], edgecolor='white', width=barWidth, label = species[5])
        # Create seven Bars
        plt.bar(r, sevenBars, bottom=[i+j+k+l+m+n for i,j,k,l,m,n in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars)], edgecolor='white', width=barWidth, label = species[6])
        # Create eight Bars
        plt.bar(r, eightBars, bottom=[i+j+k+l+m+n+o for i,j,k,l,m,n,o in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars)], edgecolor='white', width=barWidth, label = species[7])
        

        # Custom x axis
        plt.xlabel("Time (h)")

        plt.ylabel('Relative abundance (%)')
        plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)

        # Show graphic
        plt.savefig(location_of_call + '/'+ output_prefix + '-Relative_Abundance_growth.pdf', bbox_inches='tight')
        plt.close()
        plt.figure().clear()
        

    if consortia_size_wanted == 9:
        # From raw value to percentage
        totals = [i+j+k+l+m+n+o+p+q for i,j,k,l,m,n,o,p,q in zip(df2[species[0]], df2[species[1]], df2[species[2]], df2[species[3]], df2[species[4]], df2[species[5]], df2[species[6]], df2[species[7]], df2[species[8]])]
        greenBars = [i / j * 100 for i,j in zip(df2[species[0]], totals)]
        orangeBars = [i / j * 100 for i,j in zip(df2[species[1]], totals)]
        threeBars = [i / j * 100 for i,j in zip(df2[species[2]], totals)]
        fourBars = [i / j * 100 for i,j in zip(df2[species[3]], totals)]
        fiveBars = [i / j * 100 for i,j in zip(df2[species[4]], totals)]
        sixBars = [i / j * 100 for i,j in zip(df2[species[5]], totals)]
        sevenBars = [i / j * 100 for i,j in zip(df2[species[6]], totals)]
        eightBars = [i / j * 100 for i,j in zip(df2[species[7]], totals)]
        nineBars = [i / j * 100 for i,j in zip(df2[species[8]], totals)]

        # plot
        barWidth = 0.85
        # Create green Bars
        plt.bar(r, greenBars, edgecolor='white', width=barWidth, label = species[0])
        # Create orange Bars
        plt.bar(r, orangeBars, bottom=greenBars, edgecolor='white', width=barWidth, label = species[1])
        # Create three Bars
        plt.bar(r, threeBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], edgecolor='white', width=barWidth, label = species[2])
        # Create four Bars
        plt.bar(r, fourBars, bottom=[i+j+k for i,j,k in zip(greenBars, orangeBars,threeBars)], edgecolor='white', width=barWidth, label = species[3])
        # Create five Bars
        plt.bar(r, fiveBars, bottom=[i+j+k+l for i,j,k,l in zip(greenBars, orangeBars,threeBars,fourBars)], edgecolor='white', width=barWidth, label = species[4])
        # Create six Bars
        plt.bar(r, sixBars, bottom=[i+j+k+l+m for i,j,k,l,m in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars)], edgecolor='white', width=barWidth, label = species[5])
        # Create seven Bars
        plt.bar(r, sevenBars, bottom=[i+j+k+l+m+n for i,j,k,l,m,n in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars)], edgecolor='white', width=barWidth, label = species[6])
        # Create eight Bars
        plt.bar(r, eightBars, bottom=[i+j+k+l+m+n+o for i,j,k,l,m,n,o in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars)], edgecolor='white', width=barWidth, label = species[7])
        # Create nine Bars
        plt.bar(r, nineBars, bottom=[i+j+k+l+m+n+o+p for i,j,k,l,m,n,o,p in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars)], edgecolor='white', width=barWidth, label = species[8])
        

        # Custom x axis
        plt.xlabel("Time (h)")

        plt.ylabel('Relative abundance (%)')
        plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)

        # Show graphic
        plt.savefig(location_of_call + '/'+ output_prefix + '-Relative_Abundance_growth.pdf', bbox_inches='tight')
        plt.close()
        plt.figure().clear()


        
        
        
    if consortia_size_wanted == 10:
        # From raw value to percentage
        totals = [i+j+k+l+m+n+o+p+q+r for i,j,k,l,m,n,o,p,q,r in zip(df2[species[0]], df2[species[1]], df2[species[2]], df2[species[3]], df2[species[4]], df2[species[5]], df2[species[6]], df2[species[7]], df2[species[8]], df2[species[9]])]
        greenBars = [i / j * 100 for i,j in zip(df2[species[0]], totals)]
        orangeBars = [i / j * 100 for i,j in zip(df2[species[1]], totals)]
        threeBars = [i / j * 100 for i,j in zip(df2[species[2]], totals)]
        fourBars = [i / j * 100 for i,j in zip(df2[species[3]], totals)]
        fiveBars = [i / j * 100 for i,j in zip(df2[species[4]], totals)]
        sixBars = [i / j * 100 for i,j in zip(df2[species[5]], totals)]
        sevenBars = [i / j * 100 for i,j in zip(df2[species[6]], totals)]
        eightBars = [i / j * 100 for i,j in zip(df2[species[7]], totals)]
        nineBars = [i / j * 100 for i,j in zip(df2[species[8]], totals)]
        tenthBars = [i / j * 100 for i,j in zip(df2[species[9]], totals)]

        # plot
        barWidth = 0.85
        # Create green Bars
        plt.bar(r, greenBars, edgecolor='white', width=barWidth, label = species[0])
        # Create orange Bars
        plt.bar(r, orangeBars, bottom=greenBars, edgecolor='white', width=barWidth, label = species[1])
        # Create three Bars
        plt.bar(r, threeBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], edgecolor='white', width=barWidth, label = species[2])
        # Create four Bars
        plt.bar(r, fourBars, bottom=[i+j+k for i,j,k in zip(greenBars, orangeBars,threeBars)], edgecolor='white', width=barWidth, label = species[3])
        # Create five Bars
        plt.bar(r, fiveBars, bottom=[i+j+k+l for i,j,k,l in zip(greenBars, orangeBars,threeBars,fourBars)], edgecolor='white', width=barWidth, label = species[4])
        # Create six Bars
        plt.bar(r, sixBars, bottom=[i+j+k+l+m for i,j,k,l,m in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars)], edgecolor='white', width=barWidth, label = species[5])
        # Create seven Bars
        plt.bar(r, sevenBars, bottom=[i+j+k+l+m+n for i,j,k,l,m,n in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars)], edgecolor='white', width=barWidth, label = species[6])
        # Create eight Bars
        plt.bar(r, eightBars, bottom=[i+j+k+l+m+n+o for i,j,k,l,m,n,o in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars)], edgecolor='white', width=barWidth, label = species[7])
        # Create nine Bars
        plt.bar(r, nineBars, bottom=[i+j+k+l+m+n+o+p for i,j,k,l,m,n,o,p in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars)], edgecolor='white', width=barWidth, label = species[8])
        # Create tenth Bars
        plt.bar(r, tenthBars, bottom=[i+j+k+l+m+n+o+p+q for i,j,k,l,m,n,o,p,q in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars)], edgecolor='white', width=barWidth, label = species[9])
        

        # Custom x axis
        plt.xlabel("Time (h)")

        plt.ylabel('Relative abundance (%)')
        plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)

        # Show graphic
        plt.savefig(location_of_call + '/'+ output_prefix + '-Relative_Abundance_growth.pdf', bbox_inches='tight')
        plt.close()
        plt.figure().clear()
        
        
    if consortia_size_wanted == 11:
        # From raw value to percentage
        totals = [i+j+k+l+m+n+o+p+q+r+s for i,j,k,l,m,n,o,p,q,r,s in zip(df2[species[0]], df2[species[1]], df2[species[2]], df2[species[3]], df2[species[4]], df2[species[5]], df2[species[6]], df2[species[7]], df2[species[8]], df2[species[9]], df2[species[10]])]
        greenBars = [i / j * 100 for i,j in zip(df2[species[0]], totals)]
        orangeBars = [i / j * 100 for i,j in zip(df2[species[1]], totals)]
        threeBars = [i / j * 100 for i,j in zip(df2[species[2]], totals)]
        fourBars = [i / j * 100 for i,j in zip(df2[species[3]], totals)]
        fiveBars = [i / j * 100 for i,j in zip(df2[species[4]], totals)]
        sixBars = [i / j * 100 for i,j in zip(df2[species[5]], totals)]
        sevenBars = [i / j * 100 for i,j in zip(df2[species[6]], totals)]
        eightBars = [i / j * 100 for i,j in zip(df2[species[7]], totals)]
        nineBars = [i / j * 100 for i,j in zip(df2[species[8]], totals)]
        tenthBars = [i / j * 100 for i,j in zip(df2[species[9]], totals)]
        elevenBars = [i / j * 100 for i,j in zip(df2[species[10]], totals)]

        # plot
        barWidth = 0.85
        # Create green Bars
        plt.bar(r, greenBars, edgecolor='white', width=barWidth, label = species[0])
        # Create orange Bars
        plt.bar(r, orangeBars, bottom=greenBars, edgecolor='white', width=barWidth, label = species[1])
        # Create three Bars
        plt.bar(r, threeBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], edgecolor='white', width=barWidth, label = species[2])
        # Create four Bars
        plt.bar(r, fourBars, bottom=[i+j+k for i,j,k in zip(greenBars, orangeBars,threeBars)], edgecolor='white', width=barWidth, label = species[3])
        # Create five Bars
        plt.bar(r, fiveBars, bottom=[i+j+k+l for i,j,k,l in zip(greenBars, orangeBars,threeBars,fourBars)], edgecolor='white', width=barWidth, label = species[4])
        # Create six Bars
        plt.bar(r, sixBars, bottom=[i+j+k+l+m for i,j,k,l,m in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars)], edgecolor='white', width=barWidth, label = species[5])
        # Create seven Bars
        plt.bar(r, sevenBars, bottom=[i+j+k+l+m+n for i,j,k,l,m,n in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars)], edgecolor='white', width=barWidth, label = species[6])
        # Create eight Bars
        plt.bar(r, eightBars, bottom=[i+j+k+l+m+n+o for i,j,k,l,m,n,o in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars)], edgecolor='white', width=barWidth, label = species[7])
        # Create nine Bars
        plt.bar(r, nineBars, bottom=[i+j+k+l+m+n+o+p for i,j,k,l,m,n,o,p in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars)], edgecolor='white', width=barWidth, label = species[8])
        # Create tenth Bars
        plt.bar(r, tenthBars, bottom=[i+j+k+l+m+n+o+p+q for i,j,k,l,m,n,o,p,q in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars)], edgecolor='white', width=barWidth, label = species[9])
        # Create eleventh Bars
        plt.bar(r, elevenBars, bottom=[i+j+k+l+m+n+o+p+q+r for i,j,k,l,m,n,o,p,q,r in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars,tenthBars)], edgecolor='white', width=barWidth, label = species[10])
        

        # Custom x axis
        plt.xlabel("Time (h)")

        plt.ylabel('Relative abundance (%)')
        plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)

        # Show graphic
        plt.savefig(location_of_call + '/'+ output_prefix + '-Relative_Abundance_growth.pdf', bbox_inches='tight')
        plt.close()
        plt.figure().clear()
        

    if consortia_size_wanted == 12:
        # From raw value to percentage
        totals = [i+j+k+l+m+n+o+p+q+r+s+t for i,j,k,l,m,n,o,p,q,r,s,t in zip(df2[species[0]], df2[species[1]], df2[species[2]], df2[species[3]], df2[species[4]], df2[species[5]], df2[species[6]], df2[species[7]], df2[species[8]], df2[species[9]], df2[species[10]], df2[species[11]])]
        greenBars = [i / j * 100 for i,j in zip(df2[species[0]], totals)]
        orangeBars = [i / j * 100 for i,j in zip(df2[species[1]], totals)]
        threeBars = [i / j * 100 for i,j in zip(df2[species[2]], totals)]
        fourBars = [i / j * 100 for i,j in zip(df2[species[3]], totals)]
        fiveBars = [i / j * 100 for i,j in zip(df2[species[4]], totals)]
        sixBars = [i / j * 100 for i,j in zip(df2[species[5]], totals)]
        sevenBars = [i / j * 100 for i,j in zip(df2[species[6]], totals)]
        eightBars = [i / j * 100 for i,j in zip(df2[species[7]], totals)]
        nineBars = [i / j * 100 for i,j in zip(df2[species[8]], totals)]
        tenthBars = [i / j * 100 for i,j in zip(df2[species[9]], totals)]
        elevenBars = [i / j * 100 for i,j in zip(df2[species[10]], totals)]
        twelthBars = [i / j * 100 for i,j in zip(df2[species[11]], totals)]

        # plot
        barWidth = 0.85
        # Create green Bars
        plt.bar(r, greenBars, edgecolor='white', width=barWidth, label = species[0])
        # Create orange Bars
        plt.bar(r, orangeBars, bottom=greenBars, edgecolor='white', width=barWidth, label = species[1])
        # Create three Bars
        plt.bar(r, threeBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], edgecolor='white', width=barWidth, label = species[2])
        # Create four Bars
        plt.bar(r, fourBars, bottom=[i+j+k for i,j,k in zip(greenBars, orangeBars,threeBars)], edgecolor='white', width=barWidth, label = species[3])
        # Create five Bars
        plt.bar(r, fiveBars, bottom=[i+j+k+l for i,j,k,l in zip(greenBars, orangeBars,threeBars,fourBars)], edgecolor='white', width=barWidth, label = species[4])
        # Create six Bars
        plt.bar(r, sixBars, bottom=[i+j+k+l+m for i,j,k,l,m in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars)], edgecolor='white', width=barWidth, label = species[5])
        # Create seven Bars
        plt.bar(r, sevenBars, bottom=[i+j+k+l+m+n for i,j,k,l,m,n in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars)], edgecolor='white', width=barWidth, label = species[6])
        # Create eight Bars
        plt.bar(r, eightBars, bottom=[i+j+k+l+m+n+o for i,j,k,l,m,n,o in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars)], edgecolor='white', width=barWidth, label = species[7])
        # Create nine Bars
        plt.bar(r, nineBars, bottom=[i+j+k+l+m+n+o+p for i,j,k,l,m,n,o,p in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars)], edgecolor='white', width=barWidth, label = species[8])
        # Create tenth Bars
        plt.bar(r, tenthBars, bottom=[i+j+k+l+m+n+o+p+q for i,j,k,l,m,n,o,p,q in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars)], edgecolor='white', width=barWidth, label = species[9])
        # Create eleventh Bars
        plt.bar(r, elevenBars, bottom=[i+j+k+l+m+n+o+p+q+r for i,j,k,l,m,n,o,p,q,r in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars,tenthBars)], edgecolor='white', width=barWidth, label = species[10])
        # Create twleth Bars
        plt.bar(r, twelthBars, bottom=[i+j+k+l+m+n+o+p+q+r+s for i,j,k,l,m,n,o,p,q,r,s in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars,tenthBars,elevenBars)], edgecolor='white', width=barWidth, label = species[11])
        

        # Custom x axis
        plt.xlabel("Time (h)")

        plt.ylabel('Relative abundance (%)')
        plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)
        

        # Show graphic
        plt.savefig(location_of_call + '/'+ output_prefix + '-Relative_Abundance_growth.pdf', bbox_inches='tight')
        plt.close()
        plt.figure().clear()
        
        
        
        
        
    
        
    if consortia_size_wanted == 13:
        # From raw value to percentage
        totals = [i+j+k+l+m+n+o+p+q+r+s+t+x for i,j,k,l,m,n,o,p,q,r,s,t,x in zip(df2[species[0]], df2[species[1]], df2[species[2]], df2[species[3]], df2[species[4]], df2[species[5]], df2[species[6]], df2[species[7]], df2[species[8]], df2[species[9]], df2[species[10]], df2[species[11]], df2[species[12]])]
        greenBars = [i / j * 100 for i,j in zip(df2[species[0]], totals)]
        orangeBars = [i / j * 100 for i,j in zip(df2[species[1]], totals)]
        threeBars = [i / j * 100 for i,j in zip(df2[species[2]], totals)]
        fourBars = [i / j * 100 for i,j in zip(df2[species[3]], totals)]
        fiveBars = [i / j * 100 for i,j in zip(df2[species[4]], totals)]
        sixBars = [i / j * 100 for i,j in zip(df2[species[5]], totals)]
        sevenBars = [i / j * 100 for i,j in zip(df2[species[6]], totals)]
        eightBars = [i / j * 100 for i,j in zip(df2[species[7]], totals)]
        nineBars = [i / j * 100 for i,j in zip(df2[species[8]], totals)]
        tenthBars = [i / j * 100 for i,j in zip(df2[species[9]], totals)]
        elevenBars = [i / j * 100 for i,j in zip(df2[species[10]], totals)]
        twelthBars = [i / j * 100 for i,j in zip(df2[species[11]], totals)]
        thirtenthBars = [i / j * 100 for i,j in zip(df2[species[12]], totals)]

        # plot
        barWidth = 0.85
        # Create green Bars
        plt.bar(r, greenBars, edgecolor='white', width=barWidth, label = species[0])
        # Create orange Bars
        plt.bar(r, orangeBars, bottom=greenBars, edgecolor='white', width=barWidth, label = species[1])
        # Create three Bars
        plt.bar(r, threeBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], edgecolor='white', width=barWidth, label = species[2])
        # Create four Bars
        plt.bar(r, fourBars, bottom=[i+j+k for i,j,k in zip(greenBars, orangeBars,threeBars)], edgecolor='white', width=barWidth, label = species[3])
        # Create five Bars
        plt.bar(r, fiveBars, bottom=[i+j+k+l for i,j,k,l in zip(greenBars, orangeBars,threeBars,fourBars)], edgecolor='white', width=barWidth, label = species[4])
        # Create six Bars
        plt.bar(r, sixBars, bottom=[i+j+k+l+m for i,j,k,l,m in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars)], edgecolor='white', width=barWidth, label = species[5])
        # Create seven Bars
        plt.bar(r, sevenBars, bottom=[i+j+k+l+m+n for i,j,k,l,m,n in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars)], edgecolor='white', width=barWidth, label = species[6])
        # Create eight Bars
        plt.bar(r, eightBars, bottom=[i+j+k+l+m+n+o for i,j,k,l,m,n,o in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars)], edgecolor='white', width=barWidth, label = species[7])
        # Create nine Bars
        plt.bar(r, nineBars, bottom=[i+j+k+l+m+n+o+p for i,j,k,l,m,n,o,p in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars)], edgecolor='white', width=barWidth, label = species[8])
        # Create tenth Bars
        plt.bar(r, tenthBars, bottom=[i+j+k+l+m+n+o+p+q for i,j,k,l,m,n,o,p,q in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars)], edgecolor='white', width=barWidth, label = species[9])
        # Create eleventh Bars
        plt.bar(r, elevenBars, bottom=[i+j+k+l+m+n+o+p+q+r for i,j,k,l,m,n,o,p,q,r in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars,tenthBars)], edgecolor='white', width=barWidth, label = species[10])
        # Create twleth Bars
        plt.bar(r, twelthBars, bottom=[i+j+k+l+m+n+o+p+q+r+s for i,j,k,l,m,n,o,p,q,r,s in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars,tenthBars,elevenBars)], edgecolor='white', width=barWidth, label = species[11]) 
        # Create 13 Bars
        plt.bar(r, thirtenthBars, bottom=[i+j+k+l+m+n+o+p+q+r+s+t for i,j,k,l,m,n,o,p,q,r,s,t in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars,tenthBars,elevenBars,twelthBars)], edgecolor='white', width=barWidth, label = species[12])
        

        # Custom x axis
        plt.xlabel("Time (h)")

        plt.ylabel('Relative abundance (%)')
        plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)

        # Show graphic
        plt.savefig(location_of_call + '/'+ output_prefix + '-Relative_Abundance_growth.pdf', bbox_inches='tight')
        plt.close()
        plt.figure().clear()
        
        
        
        
        
        
    if consortia_size_wanted == 14:
        # From raw value to percentage
        totals = [i+j+k+l+m+n+o+p+q+r+s+t+x+y for i,j,k,l,m,n,o,p,q,r,s,t,x,y in zip(df2[species[0]], df2[species[1]], df2[species[2]], df2[species[3]], df2[species[4]], df2[species[5]], df2[species[6]], df2[species[7]], df2[species[8]], df2[species[9]], df2[species[10]], df2[species[11]], df2[species[12]], df2[species[13]])]
        greenBars = [i / j * 100 for i,j in zip(df2[species[0]], totals)]
        orangeBars = [i / j * 100 for i,j in zip(df2[species[1]], totals)]
        threeBars = [i / j * 100 for i,j in zip(df2[species[2]], totals)]
        fourBars = [i / j * 100 for i,j in zip(df2[species[3]], totals)]
        fiveBars = [i / j * 100 for i,j in zip(df2[species[4]], totals)]
        sixBars = [i / j * 100 for i,j in zip(df2[species[5]], totals)]
        sevenBars = [i / j * 100 for i,j in zip(df2[species[6]], totals)]
        eightBars = [i / j * 100 for i,j in zip(df2[species[7]], totals)]
        nineBars = [i / j * 100 for i,j in zip(df2[species[8]], totals)]
        tenthBars = [i / j * 100 for i,j in zip(df2[species[9]], totals)]
        elevenBars = [i / j * 100 for i,j in zip(df2[species[10]], totals)]
        twelthBars = [i / j * 100 for i,j in zip(df2[species[11]], totals)]
        thirtenthBars = [i / j * 100 for i,j in zip(df2[species[12]], totals)]
        forteenthBars = [i / j * 100 for i,j in zip(df2[species[13]], totals)]

        # plot
        barWidth = 0.85
        # Create green Bars
        plt.bar(r, greenBars, edgecolor='white', width=barWidth, label = species[0])
        # Create orange Bars
        plt.bar(r, orangeBars, bottom=greenBars, edgecolor='white', width=barWidth, label = species[1])
        # Create three Bars
        plt.bar(r, threeBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], edgecolor='white', width=barWidth, label = species[2])
        # Create four Bars
        plt.bar(r, fourBars, bottom=[i+j+k for i,j,k in zip(greenBars, orangeBars,threeBars)], edgecolor='white', width=barWidth, label = species[3])
        # Create five Bars
        plt.bar(r, fiveBars, bottom=[i+j+k+l for i,j,k,l in zip(greenBars, orangeBars,threeBars,fourBars)], edgecolor='white', width=barWidth, label = species[4])
        # Create six Bars
        plt.bar(r, sixBars, bottom=[i+j+k+l+m for i,j,k,l,m in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars)], edgecolor='white', width=barWidth, label = species[5])
        # Create seven Bars
        plt.bar(r, sevenBars, bottom=[i+j+k+l+m+n for i,j,k,l,m,n in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars)], edgecolor='white', width=barWidth, label = species[6])
        # Create eight Bars
        plt.bar(r, eightBars, bottom=[i+j+k+l+m+n+o for i,j,k,l,m,n,o in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars)], edgecolor='white', width=barWidth, label = species[7])
        # Create nine Bars
        plt.bar(r, nineBars, bottom=[i+j+k+l+m+n+o+p for i,j,k,l,m,n,o,p in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars)], edgecolor='white', width=barWidth, label = species[8])
        # Create tenth Bars
        plt.bar(r, tenthBars, bottom=[i+j+k+l+m+n+o+p+q for i,j,k,l,m,n,o,p,q in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars)], edgecolor='white', width=barWidth, label = species[9])
        # Create eleventh Bars
        plt.bar(r, elevenBars, bottom=[i+j+k+l+m+n+o+p+q+r for i,j,k,l,m,n,o,p,q,r in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars,tenthBars)], edgecolor='white', width=barWidth, label = species[10])
        # Create twleth Bars
        plt.bar(r, twelthBars, bottom=[i+j+k+l+m+n+o+p+q+r+s for i,j,k,l,m,n,o,p,q,r,s in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars,tenthBars,elevenBars)], edgecolor='white', width=barWidth, label = species[11]) 
        # Create 13 Bars
        plt.bar(r, thirtenthBars, bottom=[i+j+k+l+m+n+o+p+q+r+s+t for i,j,k,l,m,n,o,p,q,r,s,t in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars,tenthBars,elevenBars,twelthBars)], edgecolor='white', width=barWidth, label = species[12])   
        # Create 14 Bars
        plt.bar(r, forteenthBars, bottom=[i+j+k+l+m+n+o+p+q+r+s+t+x for i,j,k,l,m,n,o,p,q,r,s,t,x in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars,tenthBars,elevenBars,twelthBars,thirtenthBars)], edgecolor='white', width=barWidth, label = species[13])
        

        # Custom x axis
        plt.xlabel("Time (h)")

        plt.ylabel('Relative abundance (%)')
        plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)

        # Show graphic
        plt.savefig(location_of_call + '/'+ output_prefix + '-Relative_Abundance_growth.pdf', bbox_inches='tight')
        plt.close()
        plt.figure().clear()
        

            
            
            
        
    if consortia_size_wanted == 15:
        # From raw value to percentage
        totals = [i+j+k+l+m+n+o+p+q+r+s+t+x+y+z for i,j,k,l,m,n,o,p,q,r,s,t,x,y,z in zip(df2[species[0]], df2[species[1]], df2[species[2]], df2[species[3]], df2[species[4]], df2[species[5]], df2[species[6]], df2[species[7]], df2[species[8]], df2[species[9]], df2[species[10]], df2[species[11]], df2[species[12]], df2[species[13]], df2[species[14]])]
        greenBars = [i / j * 100 for i,j in zip(df2[species[0]], totals)]
        orangeBars = [i / j * 100 for i,j in zip(df2[species[1]], totals)]
        threeBars = [i / j * 100 for i,j in zip(df2[species[2]], totals)]
        fourBars = [i / j * 100 for i,j in zip(df2[species[3]], totals)]
        fiveBars = [i / j * 100 for i,j in zip(df2[species[4]], totals)]
        sixBars = [i / j * 100 for i,j in zip(df2[species[5]], totals)]
        sevenBars = [i / j * 100 for i,j in zip(df2[species[6]], totals)]
        eightBars = [i / j * 100 for i,j in zip(df2[species[7]], totals)]
        nineBars = [i / j * 100 for i,j in zip(df2[species[8]], totals)]
        tenthBars = [i / j * 100 for i,j in zip(df2[species[9]], totals)]
        elevenBars = [i / j * 100 for i,j in zip(df2[species[10]], totals)]
        twelthBars = [i / j * 100 for i,j in zip(df2[species[11]], totals)]
        thirtenthBars = [i / j * 100 for i,j in zip(df2[species[12]], totals)]
        forteenthBars = [i / j * 100 for i,j in zip(df2[species[13]], totals)]
        fithteenthBars = [i / j * 100 for i,j in zip(df2[species[14]], totals)]

        # plot
        barWidth = 0.85
        # Create green Bars
        plt.bar(r, greenBars, edgecolor='white', width=barWidth, label = species[0])
        # Create orange Bars
        plt.bar(r, orangeBars, bottom=greenBars, edgecolor='white', width=barWidth, label = species[1])
        # Create three Bars
        plt.bar(r, threeBars, bottom=[i+j for i,j in zip(greenBars, orangeBars)], edgecolor='white', width=barWidth, label = species[2])
        # Create four Bars
        plt.bar(r, fourBars, bottom=[i+j+k for i,j,k in zip(greenBars, orangeBars,threeBars)], edgecolor='white', width=barWidth, label = species[3])
        # Create five Bars
        plt.bar(r, fiveBars, bottom=[i+j+k+l for i,j,k,l in zip(greenBars, orangeBars,threeBars,fourBars)], edgecolor='white', width=barWidth, label = species[4])
        # Create six Bars
        plt.bar(r, sixBars, bottom=[i+j+k+l+m for i,j,k,l,m in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars)], edgecolor='white', width=barWidth, label = species[5])
        # Create seven Bars
        plt.bar(r, sevenBars, bottom=[i+j+k+l+m+n for i,j,k,l,m,n in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars)], edgecolor='white', width=barWidth, label = species[6])
        # Create eight Bars
        plt.bar(r, eightBars, bottom=[i+j+k+l+m+n+o for i,j,k,l,m,n,o in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars)], edgecolor='white', width=barWidth, label = species[7])
        # Create nine Bars
        plt.bar(r, nineBars, bottom=[i+j+k+l+m+n+o+p for i,j,k,l,m,n,o,p in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars)], edgecolor='white', width=barWidth, label = species[8])
        # Create tenth Bars
        plt.bar(r, tenthBars, bottom=[i+j+k+l+m+n+o+p+q for i,j,k,l,m,n,o,p,q in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars)], edgecolor='white', width=barWidth, label = species[9])
        # Create eleventh Bars
        plt.bar(r, elevenBars, bottom=[i+j+k+l+m+n+o+p+q+r for i,j,k,l,m,n,o,p,q,r in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars,tenthBars)], edgecolor='white', width=barWidth, label = species[10])
        # Create twleth Bars
        plt.bar(r, twelthBars, bottom=[i+j+k+l+m+n+o+p+q+r+s for i,j,k,l,m,n,o,p,q,r,s in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars,tenthBars,elevenBars)], edgecolor='white', width=barWidth, label = species[11]) 
        # Create 13 Bars
        plt.bar(r, thirtenthBars, bottom=[i+j+k+l+m+n+o+p+q+r+s+t for i,j,k,l,m,n,o,p,q,r,s,t in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars,tenthBars,elevenBars,twelthBars)], edgecolor='white', width=barWidth, label = species[12])   
        # Create 14 Bars
        plt.bar(r, forteenthBars, bottom=[i+j+k+l+m+n+o+p+q+r+s+t+x for i,j,k,l,m,n,o,p,q,r,s,t,x in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars,tenthBars,elevenBars,twelthBars,thirtenthBars)], edgecolor='white', width=barWidth, label = species[13])
        # Create 15 Bars
        plt.bar(r, fithteenthBars, bottom=[i+j+k+l+m+n+o+p+q+r+s+t+x+y for i,j,k,l,m,n,o,p,q,r,s,t,x,y in zip(greenBars, orangeBars,threeBars,fourBars,fiveBars,sixBars,sevenBars,eightBars,nineBars,tenthBars,elevenBars,twelthBars,thirtenthBars,forteenthBars)], edgecolor='white', width=barWidth, label = species[14])
        

        # Custom x axis
        plt.xlabel("Time (h)")

        plt.ylabel('Relative abundance (%)')
        plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)

        # Show graphic
        plt.savefig(location_of_call + '/'+ output_prefix + '-Relative_Abundance_growth.pdf', bbox_inches='tight')
        plt.close()
        plt.figure().clear()


## end of plots from step 4

raw_comparison_values = {} # Bacteria; [number of times selected, [list of positions when selected]]

for column in samples[samples_to_be_studied].columns: # For each sample        ###########
    consortia_matches, consortia_mismatches, cul_accounted,cul_mismatches, consortia = Stored_first_round_consortia[column]
    for position,sp in enumerate(consortia):
        if sp in raw_comparison_values:
            previous_values = raw_comparison_values[sp]
            previous_values[0] +=1
            previous_values[1].append(position + 1)
        else:
            raw_comparison_values[sp] = [1,[position + 1]]

            
            
genomes_coverages = {} # sp; 98.7
for genome in list(raw_comparison_values.keys()):
    matching_list = []
    for column in list(samples[samples_to_be_studied].columns): # For each sample        ###########
        accouted_pfams_percentage = 0.0
        samples_pfams = list(sample_data.loc[sample_data == 1].index)

        genome_data = genomes[genome] # This section creates the single genome file
        redundant_functional_list = list(genome_data.loc[genome_data == 1].index)
        genome_pfams = list(set(redundant_functional_list))

        match = len(list(set(samples_pfams).intersection(genome_pfams))) # Identifies the common
        accouted_pfams_percentage = (match/ float(len(samples_pfams)))*100
        matching_list.append(accouted_pfams_percentage)
    genomes_coverages[genome] = mean(matching_list)
    
    
converted_comparison_values = {}

for sp, information in raw_comparison_values.items():
    if sp in consortia_selected:
        converted_comparison_values[sp] = [(information[0]/len(samples[samples_to_be_studied].columns))*100, mean(information[1]),'Included', genomes_coverages[sp]]
    else:
        converted_comparison_values[sp] = [(information[0]/len(samples[samples_to_be_studied].columns))*100, mean(information[1]),'Not selected',genomes_coverages[sp]]
    
    
plotting_data = pd.DataFrame.from_dict(converted_comparison_values, orient='index', columns=['Selection prevalence (n = ' + str(len(samples[samples_to_be_studied].columns)) + ')', 'Average selection position','Final consortium','Average matching Pfams (%)'])





plt.figure(figsize=(8, 5))

ax = sns.scatterplot(data=plotting_data, y='Selection prevalence (n = ' + str(len(samples[samples_to_be_studied].columns)) + ')', x="Average selection position", hue="Final consortium",palette=['#14397d','#D3D3D3'], alpha=0.5, size='Average matching Pfams (%)')

ax.set(xlim=(0, 20),ylim=(0, 100))

ax.legend(loc="upper left", edgecolor="black", bbox_to_anchor=(1,1), ncol=1)


plt.savefig(location_of_call + '/'+ output_prefix + '-Consortium_selection_prevalence.pdf', bbox_inches='tight')
plt.close()
plt.figure().clear()

sorted_consortia = dict(sorted(consortia_mean_Scores.items(), key=lambda item: item[1]))

renamed_sorted = {}

num = 0
for k,v in sorted_consortia.items():
    num +=1
    if num > 0:
        renamed_sorted[num] = v
        
        
consortia_selected_key = [k for k, v in renamed_sorted.items() if v == consortia_mean_Scores[consortia_selected]]

plt.figure(figsize=(8, 5))

plt.plot(renamed_sorted.keys(), renamed_sorted.values(), 'g-',linewidth=2, color='#D3D3D3')
plt.plot(consortia_selected_key, consortia_mean_Scores[consortia_selected], 'bo',linewidth=5, color='#14397d')

plt.xlabel('Iteration')
plt.ylabel('MiMiC group-wide score')

plt.savefig(location_of_call + '/'+ output_prefix + '-Consortium_scores.pdf', bbox_inches='tight')
plt.close()
plt.figure().clear()

step4_end_time = process_time()

execution_time = (end_time - start_time) / 60 
step4_execution_time = (step4_end_time - step4_start_time) / 60 
print(f"Time taken for step 1-3: {execution_time:.2f} minutes")

if models_folder is not None and models_folder != " ":
    print(f"Time taken for step 4: {step4_execution_time:.2f} minutes")

print ('::::::::: MiMiC2 is complete.')
