#! /usr/bin/env python3

# Generates training data for an AI model that predicts the accuracy of the variant calling based on the coverage pattern.
# Generates training data by bootstrapping from the imported experimental data


import numpy as np
import pickle, os, csv, shutil
import pandas as pd


# Functions to interact with freyja   
def subsample_freyja_inputs (mask, mask_idx, file_dir='.'):    
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
        os.mkdir('%s/%d' % (file_dir, subfolder_idx))
        os.mkdir('%s/%d/depths' % (file_dir, subfolder_idx))
        os.mkdir('%s/%d/variants' % (file_dir, subfolder_idx))
        
    with open('%s/%d/depths/%d.tsv' % (file_dir, subfolder_idx, mask_idx), 'w') as outfile:
        outfile.writelines(depths2print)
       
    with open('%s/%d/variants/%d.tsv' % (file_dir, subfolder_idx, mask_idx), 'w') as outfile:
        outfile.writelines(vars2print)
    

def parse_freyja_output_old (file_name):
    freyja_raw = pd.read_table(file_name, index_col=0)
    full_var_calls = dict(eval( pd.Series(freyja_raw.loc['summarized'][0])[0].replace('inf', '0').replace('nan', '0') ))
    
    if 'Omicron' in full_var_calls.keys():
        return full_var_calls['Omicron']
    else:
        return 0


def parse_freyja_output (file_name):
    freyja_raw = pd.read_table(file_name, index_col=0)
    full_var_calls = dict(eval( pd.Series(freyja_raw.loc['lineages'][0])[0].replace('inf', '0').replace('nan', '0') ))
    
    lineages = eval( pd.Series(freyja_raw.loc['lineages'][0])[0].replace(' ', ','))
    abundances = eval( ','.join(pd.Series(freyja_raw.loc['abundances'][0])[0].split()) )

    # Calculate the total frequency of BA.1 and BA.2 sub-sublineages
    BA1_cumulated = 0
    BA2_cumulated = 0
    for (var,freq) in zip(lineage,abundances):
        if 'BA.1' in var:
            BA1_cumulated += freq
        if 'BA.2' in var:
            BA2_cumulated += freq
    
    return (BA1_cumulated, BA2_cumulated)



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


in_dir = '/hpc/scratch/Tunc.Kayikcioglu/freyja_input/'
out_dir = '/hpc/scratch/Tunc.Kayikcioglu/freyja_output/'
if os.path.exists(in_dir):
    print('Deleting the previous %s...' % in_dir)
    shutil.rmtree(in_dir)
os.mkdir(in_dir)


print('Generating subsampling mashes for the training sets...')
masks = np.empty((num_samples,num_features))
for i in range(num_samples):
    # Random number of loci at random locations are selected to be missing
    loci_present = np.ones(num_features)
    fail_ratio = np.random.rand()
    
    cursor = 0
    while cursor < num_features:
        prev_cursor = cursor
        cursor += np.random.randint(0,501)
        if np.random.rand() < fail_ratio:
            loci_present[prev_cursor:cursor] = 0
    
    subsample_freyja_inputs(loci_present, i, file_dir=in_dir)
    masks[i,:] = loci_present


print('Feeding subsampled inputs to Freyja...')
os.system("./executeFreyja.nf --inpath %s --outpath %s -w ~/scratch/work" % (in_dir, out_dir))


print('Parsing freyja output...')
freyja_truths = np.empty(num_samples)
for i in range(num_samples):
    folder_idx = np.floor(i/samples_per_batch)
    freyja_truths[i] = parse_freyja_output(file_name=out_dir+'/demix/%d' % i)


# Abundance estimates by Freyja for the full data set
full_var_call = parse_freyja_output(file_name='./inputs/%s/freyja.demix' % sample_name)


print('Exporting the training data to file...')
with open('./%s-trainingData.pkl' % sample_name, 'wb') as file:
    pickle.dump(masks, file)
    pickle.dump(freyja_truths, file)
    pickle.dump(full_var_call, file)





num_samples = 100000
samples_per_batch = 20

for sample_name in ['CFSANSMP000113119']:
    process_file(sample_name)


print('Done')

