#! /usr/bin/env python3

# Generates training data for an AI model that predicts the accuracy of the variant calling based on the coverage pattern.
# This is a form of bootstrapping from the imported experimental data

# source ../CFSANonly/prepareEnvironment.sh
# module load python/3.8.1

import numpy as np
import pickle, os, sys, csv
import pandas as pd


from getScratchPath import *
scratch_dir = getScratchPath()


# Functions to interact with freyja   
def subsample_freyja_inputs (mask, mask_idx, depth_file_present, depth_file_absent, var_file_header, var_file_contents, file_dir='.'):
    num_features = len(mask)
    depths2print = []
    for i in range(num_features):
        if mask[i]==1:
            depths2print.append(depth_file_present[i])
        else:
            depths2print.append(depth_file_absent[i])
    
    vars2print = [var_file_header]
    for i in np.where(mask==1)[0]:
        vars2print.append(var_file_contents[i])
        
    subfolder_idx = np.floor(mask_idx/samples_per_batch)
    if not os.path.exists('%s/%d' % (file_dir, subfolder_idx)):
        os.makedirs('%s/%d/depths' % (file_dir, subfolder_idx))
        os.makedirs('%s/%d/variants' % (file_dir, subfolder_idx))
        
    with open('%s/%d/depths/%d.tsv' % (file_dir, subfolder_idx, mask_idx), 'w') as outfile:
        outfile.writelines(depths2print)
       
    with open('%s/%d/variants/%d.tsv' % (file_dir, subfolder_idx, mask_idx), 'w') as outfile:
        outfile.writelines(vars2print)
    

# Parses the line 3: lineages & line 4: abundances in freyja.demix output.
# Returns the frequency of lineages that are desired to be quantified with a high accuracy.
sys.path.append('../')
from getDisplayName import import_freyja_demix
lineages_of_interest = ['BA.1', 'BA.2', 'BA.4', 'BA.5']

def parse_freyja_output (file_name):
    # Import the Freyja *.demix file
    (lineages, abundances, freyja_names) = import_freyja_demix(file_name)

    # Calculate the total frequency of BA.1 and BA.2 sub-sublineages
    cumulated_freq = np.zeros(len(lineages_of_interest))
    for (var,freq) in zip(lineages,abundances):
        for lineage in lineages_of_interest:
            if lineage in var:
                idx = lineages_of_interest.index(lineage)
                cumulated_freq[idx] += freq
                break
    
    return cumulated_freq



def import_input_files (sample_name):
    # Import the tsv and .depths file used for this computation
    print('Importing original data...')
    with open('./inputs/%s/freyja.depths' % sample_name, 'r') as infile:
        reader = csv.reader(infile, delimiter="\t")
        depth_file_present = []
        depth_file_absent = []
        for row in reader:
            depth_file_present.append('\t'.join(row) + '\n')
            row[3] = '0'
            depth_file_absent.append('\t'.join(row) + '\n')
     
    num_features = len(depth_file_present) # i.e., the genome size.
    with open('./inputs/%s/freyja.variants' % sample_name, 'r') as infile:
        reader = csv.reader(infile, delimiter="\t")
        var_file_header = '\t'.join(next(reader)) + '\n'
        var_file_contents = [ '' for i in range(num_features) ]
        for row in reader:
            pos_idx = int(row[1])
            var_file_contents[pos_idx] += '\t'.join(row) + '\n'
    
    return (depth_file_present, depth_file_absent, var_file_header, var_file_contents)



def make_freyja_dirs (sample_name):
    in_dir = '%s/%s-freyja-input/' % (scratch_dir, sample_name)
    # out_dir = '%s/%s-freyja-output/' % (scratch_dir, sample_name)
    if os.path.exists(in_dir):
        print('Deleting the previous %s...' % in_dir)
        os.system("rm -r %s" % in_dir)
        
    os.mkdir(in_dir)



# For each slurm job to be submitted, the number of sub-sampled input file pairs to be
# included per job. Setting to 1 means each file is processed on a separate job, which
# increases resource (de-)allocation overload. Setting to a too high value decreases
# the extent of parallelisation and hence the total lead time might increase.
samples_per_batch = 20


def process_file(sample_name, masks):
    # Loads the freyja-relevant files corresponding to the un-subsampled files.
    (depth_file_present, depth_file_absent, var_file_header, var_file_contents) = import_input_files(sample_name)
    
    # Prepare the directory structure to write freyja files into.
    make_freyja_dirs(sample_name)
    
    # Generate subsampled freyja files based on the masks and save new tsv files to be used as freyja input
    num_masks = len(masks)
    for i in range(0,num_masks):
        loci_present = masks[i]
        subsample_freyja_inputs(loci_present, i, depth_file_present, depth_file_absent, var_file_header, var_file_contents, scratch_dir+'/freyja_in')

    print('Feeding subsampled inputs to Freyja...')
    os.system("./executeFreyja.nf --inpath %s/freyja_in --outpath %s -w %s/work" % (scratch_dir, scratch_dir, scratch_dir))

    print('Parsing freyja output for subsampled files...')
    freyja_truths = np.empty((num_masks, len(lineages_of_interest)))
    for i in range(num_masks):
        folder_idx = np.floor(i/samples_per_batch)
        freyja_truths[i,:] = parse_freyja_output(file_name='%s/demix/%d' % (scratch_dir, i) )
    
    # Abundance estimates by Freyja for the full data set
    full_var_call = parse_freyja_output(file_name='./inputs/%s/freyja.demix' % sample_name)
    
    print('Exporting the training data to file...')
    with open('%s/%s_training_data.pkl' % (scratch_dir, sample_name), 'wb') as file:
        pickle.dump(freyja_truths, file)
        pickle.dump(full_var_call, file)
        


# Set of binary vectors (i.e. 0-1 matrix).
# 0: genomic coordinate was not sequenced
# 1: genomic coordinate was sequenced at least 1X
def generate_random_masks(mask_file='./masks.pkl', num_masks=100, ncores=20):
    print('Generating subsampling mashes for the training sets...')
    os.system("srun -c %d ./generateMasks.py %s %d %d" % (ncores, mask_file, num_masks, ncores))
    with open(mask_file, 'rb') as file:
        masks = pickle.load(file)
    return masks
    

masks = generate_random_masks(mask_file='%s/masks.pkl' % scratch_dir, num_masks=10000)
for sample_name in ['CFSANSMP000113119']:
    process_file(sample_name, masks)


print('Done')

