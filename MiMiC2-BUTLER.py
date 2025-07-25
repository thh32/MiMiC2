#!/usr/bin/env python3.11
import glob
import pandas as pd
import operator
import collections
import numpy as np
import sys
import scipy.stats as stats
import itertools
import subprocess
import random
from tqdm import tqdm
import pandas
from statistics import mean
import os
import argparse


location_of_file = os.path.dirname(os.path.realpath(__file__))


location_of_call = os.getcwd()

version = '2024-08-08'


# User input
parser = argparse.ArgumentParser(description='MiMiC2-BUTLER v' + version)
# Options for Sample data
parser.add_argument('-s','--samples', metavar='{INPUT}', required=True, help='Provide a folder which contains all of your Pfam annotated genomes/metagenomes.')

# Option for Function database
parser.add_argument('-p','--pfam', metavar='{INPUT}',required=True, help='Pfam file e.g. Pfam-A.clans.csv, provided for Pfam v32 in `datasets/core/`')

# State tool used  
parser.add_argument('-t','--tool', metavar='{TEXT}', help='State the tool used to annotate the geomes against the Pfam database: `hmmsearch` or `hmmscan`')

# Options for output
parser.add_argument('-o','--output', metavar='{OUTPUT}', required=True, help='Prefix for all the Pfam-profile file e.g. HuSynCom.')

# Remove file ending
parser.add_argument('-e','--extension', metavar='{TEXT}', default='.hmmer', required=True, help='Provide the extension for your Pfam annotation files.')

# Remove file ending
parser.add_argument('-m','--merge', default=False, action=argparse.BooleanOptionalAction, help='Option to merge premade vector files within a folder, define -e as "profile.txt""')


args, unknown = parser.parse_known_args()



print (': Reading in the users options.')

print (args.merge)

pfam_file = args.pfam

genome_folder = args.samples

file_ending = args.extension

output_file = args.output + '-profile.txt'



pfams = []

for line in tqdm(open(pfam_file,'r')):
    timber = line.split('\t')
    pfams.append(timber[0])
    
print (': The number of pfams studied are: ' + str(len(pfams)))


print (': Users options accepted.')

counting = 0
for cfile in tqdm(glob.glob(genome_folder + '/*' + file_ending)):
	counting +=1

print (': Files meeting criteria in folder:' + str(counting))


if args.merge == True:
	print (':: Merging your profiles.')
	combined = ''
	file_num = 0
	for cfile in tqdm(glob.glob(genome_folder + '/*' + file_ending)):
		df = pd.read_csv(cfile,sep='\t')
		file_num +=1
		if file_num == 1:
			combined = df
		else:
			combined = pd.merge(combined, df, on='PfamID')
	combined.to_csv(location_of_call + '/' + output_file, sep='\t', index=False, header=True)
	sys.exit()



if args.tool == 'hmmscan':
	print (':: Handling your files and cataloging their Pfam presence/absence.')

	samples_data = {}

	for cfile in tqdm(glob.glob(genome_folder + '/*' + file_ending)):
	    #print (cfile)
	    indiv_pfam = []
	    for line in open(cfile):
	        if line.startswith('#'):
	            lolp = 0
	        else:
	            
	            pfam = list(filter(None,line.split(' ')))[1].split('.')[0]
	            if pfam in pfams:
	                indiv_pfam.append(pfam)
	            else:
	                print ('ERROR: PFAM not found in database; ' + str(pfam))
	                print ('HINT: Are you using the same Pfam versions for each step in the analysis?')
	                print ('ACTION: Forced exiting now.')
	                sys.exit()
	                
	    samples_data[cfile.split('/')[-1:][0].replace(file_ending,'')] = list(dict.fromkeys(indiv_pfam))
	    
	    
	    
	    
	outputting = open(location_of_call + '/' + output_file,'w')


	pfam_lines = {}

	header = 'PfamID'


	for pfam in pfams: # Provide an entry for every pfam so each will be accounted for
	    pfam_lines[pfam] = pfam
	    
	for sample, data in samples_data.items():
	    theader = header + '\t' + sample # add sample name to header
	    header = theader
	    
	    tpfam_lines = {}
	    for pfam, existingline in pfam_lines.items(): # Loop over every pfam so all are accounted for
	        if pfam in data: # if the pfam is in the sample add a 1 
	            tpfam_lines[pfam] = existingline + '\t1'
	        else: # if the pfam is not in the sample add a 0
	            tpfam_lines[pfam] = existingline + '\t0'
	    pfam_lines = tpfam_lines
	    
	outputting.write(header + '\n')

	for k, v in pfam_lines.items():
	    outputting.write(v + '\n')
	    
	outputting.close()


elif args.tool == 'hmmsearch':
	print (':: Handling your files and cataloging their Pfam presence/absence.')

	samples_data = {}

	for cfile in tqdm(glob.glob(genome_folder + '/*' + file_ending)):
	    #print (cfile)
	    indiv_pfam = []
	    for line in open(cfile):
	        if line.startswith('#'):
	            lolp = 0
	        else:
	            
	            pfam = list(filter(None,line.split(' ')))[3].split('.')[0]
	            if pfam in pfams:
	                indiv_pfam.append(pfam)
	            else:
	                print ('ERROR: PFAM not found in database; ' + str(pfam))
	                print ('HINT: Are you using the same Pfam versions for each step in the analysis?')
	                print ('ACTION: Forced exiting now.')
	                sys.exit()
	                
	    samples_data[cfile.split('/')[-1:][0].replace(file_ending,'')] = list(dict.fromkeys(indiv_pfam))
	    
	    
	    
	    
	outputting = open(location_of_call + '/' + output_file,'w')


	pfam_lines = {}

	header = 'PfamID'


	for pfam in pfams: # Provide an entry for every pfam so each will be accounted for
	    pfam_lines[pfam] = pfam
	    
	for sample, data in samples_data.items():
	    theader = header + '\t' + sample # add sample name to header
	    header = theader
	    
	    tpfam_lines = {}
	    for pfam, existingline in pfam_lines.items(): # Loop over every pfam so all are accounted for
	        if pfam in data: # if the pfam is in the sample add a 1 
	            tpfam_lines[pfam] = existingline + '\t1'
	        else: # if the pfam is not in the sample add a 0
	            tpfam_lines[pfam] = existingline + '\t0'
	    pfam_lines = tpfam_lines
	    
	outputting.write(header + '\n')

	for k, v in pfam_lines.items():
	    outputting.write(v + '\n')
	    
	outputting.close()


print ('::: Your Pfam profile has been prepared for you.')

