#!/usr/bin/env python3

import pickle
import sys

if len(sys.argv) != 2:
    raise Exception('Incorrect call to the script.')


infilename = sys.argv[1]

# Import the pre-processed variant definitions from file
with open(infilename, 'rb') as file:
    uniqueVarNames = pickle.load(file)
    # These can be skipped as unneeded
    # uniqueMutationLabels = pickle.load(file)
    # var2mut = pickle.load(file)
    # mut2var = pickle.load(file) # Skipped these for efficiency
    # importantVars = pickle.load(file)
    # pos2gene = pickle.load(file)
    # gene2pos = pickle.load(file)
    # sigMutationMatrix = pickle.load(file)


# Print a concise list of variants available among the pickled definitions
print('%s' % ', '.join(uniqueVarNames))
