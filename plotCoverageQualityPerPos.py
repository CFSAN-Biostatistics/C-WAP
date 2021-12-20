#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import csv
import sys


pileupFilename = sys.argv[1]
outputDirectory = sys.argv[2]

GENOME_SIZE = 29903
quality = np.zeros(GENOME_SIZE)
readDepth = np.zeros(GENOME_SIZE)

# Import the pile-up file and record the coverage and average quality per position
with open(pileupFilename) as infile:
    # Uncomment below in case the -d <maxDepth> was unset during mpileup step.
    csv.field_size_limit(int(1e7))
    reader = csv.reader(infile, delimiter="\t", quoting=csv.QUOTE_NONE)
    for row in reader:
        pos = int(row[1])-1
        if pos >= GENOME_SIZE:
            print(
                'WARNING: Pileup file contains more rows than the genome length. Omitted extra fields.')
            break

        currentDepth = int(row[3])
        if currentDepth > 0:
            readDepth[pos] = currentDepth

            # PHRED to INT conversion rule:
            # def phred2int(x): return ord(x)-33
            quality[pos] = np.mean([ord(letter) for letter in row[5]]) - 33


# Import the list of uncovered genome regions due to kit design
gapStart = np.empty(GENOME_SIZE)
gapEnd = np.empty(GENOME_SIZE)
gapCounter = 0
gapfilename = sys.argv[3]
gapfile = open(gapfilename, 'r')
for line in gapfile.readlines():
    if len(line) > 1:
        gapStart[gapCounter], gapEnd[gapCounter] = line.split()
        gapCounter += 1
gapfile.close()


posIdx = np.arange(1, GENOME_SIZE+1, 1)
FDAblue = (0, 124/255, 186/255)  # RGB color representation of the logo


#################################################################
# Generate a plot for quality vs pos and save in a file
plt.rcParams.update({'font.size': 14})

plt.plot(posIdx, quality, '.', color=FDAblue)
plt.plot(posIdx[quality < 30], quality[quality < 30], '.', color='k')

plt.xlabel('Genome position (kb)')
plt.ylabel('Average read quality')
plt.xlim([0, GENOME_SIZE+1])
plt.xticks(np.arange(0, GENOME_SIZE, 5000), ["%d" % (
    x/1000) for x in np.arange(0, GENOME_SIZE, 5000)])


# Add shading for uncovered regions
ymax = plt.gca().get_ylim()[1]
for i in range(0, gapCounter):
    plt.fill([gapStart[i], gapEnd[i], gapEnd[i], gapStart[i]], [
             0, 0, ymax, ymax], 'r', alpha=1, edgecolor='none')
    plt.text((gapStart[i]+gapEnd[i])/2-GENOME_SIZE/60,
             1.01*ymax, '*', color='r', weight='bold', size=20)

plt.ylim([0, ymax])

plt.savefig(outputDirectory + '/quality.png', dpi=200)
plt.close()


#################################################################
# Generate a 1D histogram for the quality factors and save in a file
plt.hist(quality, bins=50, color=FDAblue, edgecolor=None, density=True)
plt.xlim([0, max(40, *quality)*1.1])

(ymin, ymax) = plt.gca().get_ylim()
plt.plot([30, 30], [0, ymax], '--', color='k')
plt.ylim(ymin, ymax)

pctBelowThreshold = 100*len([x for x in quality if x < 30]) / len(quality)
(xmin, xmax) = plt.gca().get_xlim()
if (30-xmin) < 0.15*(xmax-xmin):
    plt.text(xmin, 1.05*ymax, "%.1f%%" %
             pctBelowThreshold,  color='r', size=16)
else:
    plt.text(10, 1.05*ymax, "%.1f%%" % pctBelowThreshold,  color='r', size=16)
plt.text((30+xmax)/2-5, 1.05*ymax, "%.1f%%" %
         (100-pctBelowThreshold),  color='g', size=16)

plt.xlabel('Average read quality')
plt.ylabel('Frequency')
plt.yticks([])

plt.savefig(outputDirectory + '/qualityHistogram.png', dpi=200)
plt.close()


#################################################################
# Generate a plot for sequence coverage vs pos and save in a file
# If coverage is too high, scale the axes for better view
if np.mean(readDepth) > 500:
    plt.plot(posIdx, readDepth/1000, '.', color=FDAblue)
    plt.plot(posIdx[readDepth < 100],
             readDepth[readDepth < 100]/1000, '.', color='k')
    plt.ylabel('Coverage depth (1000)')
else:
    plt.plot(posIdx, readDepth, '.', color=FDAblue)
    plt.ylabel('Coverage depth')

plt.xlabel('Genome position (kb)')
plt.xlim([0, GENOME_SIZE+1])
plt.xticks(np.arange(0, GENOME_SIZE, 5000), ["%d" % (
    x/1000) for x in np.arange(0, GENOME_SIZE, 5000)])

# Add shading for uncovered regions
ymax = plt.gca().get_ylim()[1]
for i in range(0, gapCounter):
    plt.fill([gapStart[i], gapEnd[i], gapEnd[i], gapStart[i]], [
             0, 0, ymax, ymax], 'r', alpha=1, edgecolor='none')
    plt.text((gapStart[i]+gapEnd[i])/2-GENOME_SIZE/60,
             1.01*ymax, '*', color='r', weight='bold', size=20)

plt.ylim([0, ymax])
plt.savefig(outputDirectory + '/coverage.png', dpi=200)
plt.close()


#################################################################
# Generate a histogram for the sequencing depth and save in a file
pctBelowThreshold = 100*len([x for x in readDepth if x < 100]) / len(quality)
if np.mean(readDepth) > 1000:
    plt.hist(readDepth/1000, bins=50, color=FDAblue,
             edgecolor=None, density=True)
    (ymin, ymax) = plt.gca().get_ylim()
    plt.plot([0.1, 0.1], [0, ymax], '--', color='k')
    plt.ylim(ymin, ymax)
    plt.xlabel('Average read depth (1000)')

    (xmin, xmax) = plt.gca().get_xlim()
    if (0.1-xmin) < 0.15*(xmax-xmin):
        plt.text(xmin, 1.05*ymax, "%.1f%%" %
                 pctBelowThreshold,  color='r', size=16)
    else:
        plt.text(0.8*(xmin+0.1)/2, 1.05*ymax, "%.1f%%" %
                 pctBelowThreshold,  color='r', size=16)
    plt.text(max(0.09, 0.8*(xmax+0.1)/2), 1.05*ymax, "%.1f%%" %
             (100-pctBelowThreshold),  color='g', size=16)
else:
    plt.hist(readDepth, bins=50, color=FDAblue, edgecolor=None, density=True)
    (ymin, ymax) = plt.gca().get_ylim()
    plt.plot([100, 100], [0, ymax], '--', color='k')
    plt.ylim(ymin, ymax)
    plt.xlabel('Average read depth')

    (xmin, xmax) = plt.gca().get_xlim()
    if (100-xmin) < 0.15*(xmax-xmin):
        plt.text(xmin, 1.05*ymax, "%.1f%%" %
                 pctBelowThreshold,  color='r', size=16)
    else:
        plt.text(0.8*(xmin+100)/2, 1.05*ymax, "%.1f%%" %
                 pctBelowThreshold,  color='r', size=16)
    plt.text(max(90, 0.8*(xmax+100)/2), 1.05*ymax, "%.1f%%" %
             (100-pctBelowThreshold),  color='g', size=16)

plt.ylabel('Frequency')
plt.yticks([])
plt.savefig(outputDirectory + '/depthHistogram.png', dpi=200)
plt.close()


#################################################################
# Export a csv file for pos; coverage; quality
outfilename = outputDirectory + "/pos-coverage-quality.tsv"
with open(outfilename, 'w') as outfile:
    writer = csv.writer(outfile, delimiter="\t")
    for i in range(0, GENOME_SIZE):
        writer.writerow(["%d" % posIdx[i], "%d" %
                        readDepth[i], '%.2f' % quality[i]])
