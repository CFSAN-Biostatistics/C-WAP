#!/usr/bin/env python3

# To improve performance, the mutation signatures are pre-processed and pickled
# in this script. The pickles will then be used for repetitive actions.

import numpy as np
import json
import glob
import pickle
import os

# Pre-process the genomic features of covid genome.
# Generate a mapping between ORF/gene names and genomic coordinates.
# Each entry in the dict is name->(startPos, endPos)
# constellationsDir = "/projects/covidtrakr/software/constellations/constellations"
constellationsDir = "/projects/covidtrakr/software/miniconda3/envs/pangolin/lib/python3.8/site-packages/constellations"

file = open(constellationsDir + "/data/SARS-CoV-2.json", 'r')
fileContents = json.load(file)
file.close()

gene2pos = dict()
gene2pos["NUC"] = (1, 29903)
gene2pos["DEL"] = (1, 29903)
gene2pos["ORF1AB"] = (266, 21555)
gene2pos["1AB"] = (266, 21555)
for gene in fileContents['genes']:
    genename = gene.upper()
    start = fileContents['genes'][gene]['coordinates']['from']
    end = fileContents['genes'][gene]['coordinates']['to']
    gene2pos[genename] = (start, end)

    # Also add the name version without ORF label
    if "ORF" in genename:
        gene2pos[genename[3:]] = (start, end)


# If the protein info is not already in the dictionary, record as a separate entry.
for protein in fileContents['proteins']:
    # protname = fileContents['proteins'][protein]['name'].upper()
    protname = protein.upper()
    if protname not in gene2pos.keys():
        ORFname = fileContents['proteins'][protein]['gene'].upper()
        protStart = fileContents['proteins'][protein]['coordinates']['from']
        protEnd = fileContents['proteins'][protein]['coordinates']['to']

        (geneStart, geneEnd) = gene2pos[ORFname]
        start = geneStart + 3*protStart - 3
        end = geneStart + 3*protEnd
        gene2pos[protname] = (start, end)


# Re-structure the data such that each nucleotide position is assigned
# the relevant gene/ORF/protein names
pos2gene = dict()
for i in range(1, 29904, 1):
    pos2gene[i] = list([])

for name in gene2pos:
    (start, end) = gene2pos[name]
    for i in range(max(1, start), min(end+1, 29904)):
        pos2gene[i].append(name)


# Re-compiles the obtained covid variant signatures
# Import signature database file
mut2var = dict()
var2mut = dict()
importantVars = {"VOC": [], "VOI": [], "VUI": []}

# Source of the data to be processed:
# https://github.com/cov-lineages/constellations/tree/main/constellations/definitions
for filename in glob.iglob(constellationsDir + "/definitions/*.json"):
    print(filename)
    file = open(filename, 'r')
    fileContents = json.load(file)
    file.close()

    # Disabled WHO name substitution to make omicron subvariants visible. TK March 2022.
    #if "WHO_label" in fileContents['variant']:
    #    variantName = fileContents['variant']['WHO_label']
    #else:
    # WHO label not available
    pathEndPos = filename.rindex("/")
    variantName = filename[pathEndPos+2:-5]

    varMutList = []
    for mutation in fileContents['sites']:
        mutation = mutation.upper()

        # There is inconsistent labelling in constellations S vs SPIKE as of September 2021.
        mutation = mutation.replace('SPIKE', 'S')

        # Correct formatting of the mutation designations
        # If "ORF" is not part of "ORF8", add it back.
        newmut = ""
        if mutation[0].isdigit():
            newmut = "ORF" + mutation

        if mutation[-1] == '-':
            # This is a deletion mutation, written in this format: S:TGDGD546-
            # Convert to genomic coordinates as "DEL:position:#deletedBases"
            (ORF, point) = mutation.split(':')
            residueStartPos = int(''.join(x for x in point if x.isdigit()))
            if ORF == 'NUC':
                genomicPos = residueStartPos
                deletedResidues = ''.join(x for x in point if x.isalpha())
                numDeletedBases = len(deletedResidues)
            else:
                genomicPos = 3*residueStartPos-3 + gene2pos[ORF][0]
                deletedResidues = ''.join(x for x in point if x.isalpha())
                numDeletedBases = 3*len(deletedResidues)
            newmut = "DEL:" + str(genomicPos) + ":" + str(numDeletedBases)
        else:
            # If a name with NSP/ORF1A/ORF1B was included, convert to ORF1AB
            if "NSP" in mutation or "ORF1A:" in mutation or "ORF1B:" in mutation:
                # or just a simple residue substitution
                (ORF, point) = mutation.split(':')
                refAllele = point[0]
                newAllele = point[-1]
                pos = int(point[1:-1])
                ORFstartPos = gene2pos[ORF][0]
                ORF1ABstartPos = gene2pos["ORF1AB"][0]
                newPos = int(pos + (ORFstartPos-ORF1ABstartPos)/3)
                newmut = "ORF1AB:" + refAllele + str(newPos) + newAllele

        if len(newmut) > 0:
            print('%s is equivalent to %s' % (mutation, newmut))
            mutation = newmut

        # Record the variants in the dictionary
        varMutList.append(mutation)
        if mutation in mut2var:
            mut2var[mutation].append(variantName)
        else:
            mut2var[mutation] = [variantName]
    var2mut[variantName] = varMutList

    # If this variant has been designated as a VOC etc., record this label.
    # Note that a variant can be classified in multiple categories, but the highest one will be used
    #for tag in fileContents['tags']:
    #    if "VOC" in tag:
    #        importantVars["VOC"].append(variantName)
    #        break
    #    if "VOI" in tag:
    #        importantVars["VOI"].append(variantName)
    #        break
    #    if "VUI" in tag:
    #        importantVars["VUI"].append(variantName)
    #        break
    
    # Alternatively, manual decision based on common variants in press.
    if variantName.upper() in ['B.1.617.2', 'AY.4', 'AY.4.2', 'BA.1', 'BA.2', 'BA.3']:
        importantVars["VOC"].append(variantName)
    
# Include the wt sequence in the list, which has no mutations at all (mostly for normalisation purposes)
# wt is not a member of VOC, VOI or VUI
var2mut['wt'] = []


uniqueMutationLabels = list(mut2var.keys())
uniqueMutationLabels.sort()
uniqueVarNames = list(var2mut.keys())
uniqueVarNames.sort()
NUM_MUTATIONS = len(uniqueMutationLabels)
NUM_VARIANTS = len(uniqueVarNames)

sigMutationMatrix = np.zeros((NUM_MUTATIONS+1, NUM_VARIANTS))
for mutID in range(NUM_MUTATIONS):
    mutLabel = uniqueMutationLabels[mutID]
    for varLabel in mut2var[mutLabel]:
        varID = uniqueVarNames.index(varLabel)
        sigMutationMatrix[mutID, varID] = 1

# Tle last row of the signature matrix is filled with 1000s to serve for normalisation constraint
sigMutationMatrix[-1, :] = 1000

# Save the calculated mappings in a binary for future use.
for key in importantVars.keys():
    importantVars[key].sort()

with open('./covidRefSequences/varDefinitions.pkl', 'wb') as file:
    pickle.dump(uniqueVarNames, file)
    pickle.dump(uniqueMutationLabels, file)
    pickle.dump(var2mut, file)
    pickle.dump(mut2var, file)
    pickle.dump(importantVars, file)
    pickle.dump(pos2gene, file)
    pickle.dump(gene2pos, file)
    pickle.dump(sigMutationMatrix, file)


print('Mutation definitions were imported for:')
print(uniqueVarNames)
print(importantVars)
