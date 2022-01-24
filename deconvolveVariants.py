#!/usr/bin/env python3

import numpy as np
import pickle
import csv
import sys
from sklearn.linear_model import LinearRegression


# Variant file name is passed as a command line argument
variantFilename = sys.argv[1]
outputDirectory = sys.argv[2]
variantDBfilename = sys.argv[3]


# Import the pre-processed variant definitions from file
with open(variantDBfilename, 'rb') as file:
    uniqueVarNames = pickle.load(file)
    uniqueMutationLabels = pickle.load(file)
    var2mut = pickle.load(file)
    mut2var = pickle.load(file)  # Skipped these for efficiency
    importantVars = pickle.load(file)
    pos2gene = pickle.load(file)
    gene2pos = pickle.load(file)
    sigMutationMatrix = pickle.load(file)

# Re-compute some useful parameters
NUM_MUTATIONS = len(uniqueMutationLabels)
NUM_VARIANTS = len(uniqueVarNames)


# Parse the variant call file (vcf from samtools or tsv from iVar) and
# retain the statistically highly significant mutations detected.
outfile = open(outputDirectory+"/mutationTable.html", 'w')
with open(variantFilename, encoding='cp1252') as infile:
    reader = csv.reader(infile, delimiter="\t")
    infile.readline()  # Skip header line
    # Loop over all detected variants and tabulate their identity.
    # Only consider variants that are statistically highly significant
    freqVec = np.zeros(NUM_MUTATIONS+1)
    prevMutPos = 0
    for row in reader:
        pValue = float(row[12])
        if pValue < 0.01:
            mutationPos = int(row[1])
            # The length of the "assembled" reads is longer than the reference genome. Ignore the rest.
            if mutationPos > 29903:
                break
            
            # Avoid duplicate reporting, which happens if multiple overlapping
            # ORF's are present in the GFF file provided to iVar.
            if prevMutPos == mutationPos:
                continue
            else:
                prevMutPos = mutationPos

            # Convert genomic coordinates to CDS coordinates
            if len(pos2gene[mutationPos]) >= 3:
                # This is in a known CDS
                ORFname = pos2gene[mutationPos][2]

                # This is in ORF1B which has -1 shift due to ribosome slippage
                if ORFname == "ORF1AB" and "ORF1B" in pos2gene[mutationPos]:
                    mutationPosORF = (mutationPos-gene2pos["ORF1B"][0])//3 + 1
                else:
                    mutationPosORF = (mutationPos-gene2pos[ORFname][0])//3 + 1
            else:
                # This is in a UTR
                ORFname = "NUC"
                mutationPosORF = mutationPos

            # Parse mutation data from file
            refBase = row[2]
            altBase = row[3]
            altFreq = float(row[10])
            refAA = row[16]
            altAA = row[18]

            # Construct a mutation label based on the reporting format of
            # the mutations in constellations
            if (len(altAA) == 1) and (ORFname != "NUC") and (refAA != altAA):
                # Mutations causing residue substitutions in CDS
                mutationLabel = ''.join(
                    (ORFname, ':', refAA, str(mutationPosORF), altAA))
            else:
                if (len(altBase) == 1) or (altBase[0] == '+'):
                    # Insertions (+A +AT etc.) OR substitutions in UTR OR
                    # silent residue substitutions in CDS
                    mutationLabel = ''.join(
                        ('NUC:', refBase, str(mutationPos), altBase))
                else:
                    # Deletions in CDS or UTR
                    mutationLabel = ''.join(
                        ('DEL:', str(mutationPos), ':', str(len(altBase)-1)))

            if mutationLabel in uniqueMutationLabels:
                mutIdx = uniqueMutationLabels.index(mutationLabel)
                freqVec[mutIdx] = altFreq

            # Generate html output to be used for the "detected mutations" section
            if mutationLabel in uniqueMutationLabels:
                compatibleVars = []
                for varname in mut2var[mutationLabel]:
                    compatibleVars.append(
                        "<a href=\"https://outbreak.info/situation-reports?pango="
                        + varname + "\">" + varname + "</a>")
                compatibleVars = ', '.join(compatibleVars)
            else:
                compatibleVars = "None found"

            # If present at a significant level, then include among the mutations
            # to be shown on the report file.
            if altFreq >= 0.05:
                outfile.write("<tr>\n")
                outfile.write("    <td>%d</td> <td>%s</td> <td>%s</td> <td>%.3f</td> <td nowrap>%.2E</td>"
                    "<td>%s</td> <td>%s</td>\n"
                              % (mutationPos, refBase, altBase, altFreq, pValue, mutationLabel, compatibleVars))
                outfile.write("</tr>\n")
outfile.close()

# Tle last entry of vector is set to 1000 to serve for normalisation constraint
freqVec[-1] = 1000

# Perform linear regression to estimate variant prevalences
linReg = LinearRegression(positive=True).fit(sigMutationMatrix, freqVec)
variantPercentages = 100*linReg.coef_
R2 = linReg.score(sigMutationMatrix[0:-1, :], freqVec[0:-1])

# Log the deduced variant decomposition to file for future use.
# Format: varName %present per each line
with open(outputDirectory + "/linearDeconvolution_abundance.csv", 'w') as logfile:
    for i in range(NUM_VARIANTS):
        logfile.write('%s %.1f\n' % (uniqueVarNames[i], variantPercentages[i]))


# Print the name of the most abundant variant and its fraction to the terminal.
highestFreq = np.amax(variantPercentages)
highestIdx = np.argmax(variantPercentages)
print('%.2f %s %.2f' % (highestFreq, uniqueVarNames[highestIdx], R2))


# Check if the presence of a variant is supported or not.
# This is done by verifying if all of the mutations of the variant
# have been observed. If at least one of the characteristic mutations
# is undetected, the variant is said to be not present.

def isVarSupported(varname):
    if varname in uniqueVarNames:
        varID = uniqueVarNames.index(varname)

        relevantMuts = {'supporting': [], 'unsupporting': []}
        for i in range(NUM_MUTATIONS):
            if sigMutationMatrix[i, varID] == 1:
                if freqVec[i] < 0.001:
                    # At least one of the mutations does not exist (below threshold presence)
                    relevantMuts['unsupporting'].append(
                        uniqueMutationLabels[i])
                else:
                    # This mutation is present to sufficient degree.
                    relevantMuts['supporting'].append(uniqueMutationLabels[i])

    else:
        raise LookupError(
            "The mutation signatures for variant %s not found" % var)
    return relevantMuts


outfile = open(outputDirectory + '/VOC-VOIsupportTable.html', 'w')
outfile.write("""<p> Under the assumption that the presence of a variant requires the detection
    of all respective mutations of the variant, the characteric mutations which support the presence of
    the respective variant are indicated in the respective column of the table. Numbers show the number
    of mutations detected, if any, and the number of mutations expected to be present based on the
    variant definitions.\n""")


# Build one html table for each genre of variants (VOC/VOI/VUI)
def buildVarSupportTable(varGenre):
    outfile.write("<table>\n<tr>\n    <td>%s</td>" % varGenre)
    for varname in importantVars[varGenre]:
        outfile.write("<td>%s</td>" % ("<a href=\"https://outbreak.info/situation-reports?pango=" +
                      varname + "\">" + varname + "</a>"))

    outfile.write(
        "</tr>\n<tr>\n    <td>Characteristic mutations <br> detected</td>")
    for varname in importantVars[varGenre]:
        mutSupport = isVarSupported(varname)
        outfile.write("<td>")
        outfile.write('(%d of %d)<br>' % (len(mutSupport['supporting']),
                len(mutSupport['supporting'])+len(mutSupport['unsupporting'])))
        for mut in mutSupport['supporting']:
            outfile.write('%s<br>' % mut)
        outfile.write("</td>")
    outfile.write("\n</tr>\n</table>\n<br>\n")


if len(importantVars["VOC"]) > 0:
    buildVarSupportTable("VOC")

if len(importantVars["VOI"]) > 0:
    buildVarSupportTable("VOI")

if len(importantVars["VUI"]) > 0:
    buildVarSupportTable("VUI")

outfile.close()
