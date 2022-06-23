#! /usr/bin/env python3

# Generates random masks simulating presence/absence of coverage.
# Writes out into a pkl file as a single matrix.
# Execution is parallel across all cores, so a compute node is recommended if available.

import sys, pickle
import numpy as np
from multiprocessing import Pool, Lock, shared_memory


if len(sys.argv) != 4:
    raise Exception('Incorrect call to the script.')

mask_file = sys.argv[1]
num_masks = int(sys.argv[2])
ncores = int(sys.argv[3])
num_features = 29903 # SC2 genome length


def make_mask (mask_id):
    # Random number of loci at random locations are selected to be missing
    loci_present = np.ones(num_features)
    coverage_ratio = np.random.rand() # % genome covered, needs tuning
        
    cursor = 0
    while cursor < num_features:
        prev_cursor = cursor
        cursor += np.random.randint(0,501) # Range of gap lengths, needs tuning
        if np.random.rand() < coverage_ratio:
            loci_present[prev_cursor:cursor] = 0
    
    return loci_present

 
# Spawn multiple processes
pool = Pool(processes=ncores)
masks = pool.map(make_mask, range(num_masks))


# Wait for all processes to complete
pool.close()
pool.join()


# Record the simulated masks to a file
masks = np.array(masks)
with open(mask_file, 'wb') as file:
    pickle.dump(masks, file)

