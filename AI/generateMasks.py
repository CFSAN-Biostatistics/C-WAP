#! /usr/bin/env python3

# Generates random masks simulating presence/absence of coverage.
# Writes out into a pkl file as a single matrix.
# Execution is parallel across all cores, so a compute node is recommended if available.

import sys, pickle
import numpy as np
from multiprocessing import Pool, Lock, shared_memory



if len(sys.argv) != 3:
    raise Exception('Incorrect call to the script.')

mask_file = sys.argv[1]
num_total_masks = sys.argv[2]
num_features = 29903 # SC2 genome length


def make_mask ():
    masks = np.empty((1,num_features))
    for i in range(num_masks):
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
		
 
# Spawn multiple processes
num_processes = 20
num_masks_per_core = 
pool = Pool(processes=20, maxtasksperchild=1)
masks = pool.map(makeTree, range(num_trees))
all_masks = 


# Wait for all processes to complete
pool.close()
pool.join()


with open(mask_file, 'wb') as file:
    pickle.dump(all_masks, file)
