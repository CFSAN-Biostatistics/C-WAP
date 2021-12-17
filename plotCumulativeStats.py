#!/usr/bin/env python3
# Plots the QC metrics of all samples on consolidates 2D histograms etc. to guide an overall judgement.


import numpy as np
import matplotlib.pyplot as plt
import csv
import datetime


########################################################################

class Sample:
    def __init__(self, srr, collectionDate, collectionSite):
        self.srr = srr
        (YYYY, MM, DD) = collectionDate.split('-')
        self.collectionDate = datetime.datetime(int(YYYY), int(MM), int(DD))
        self.collectionSite = collectionSite


# Import the list of relevant SRRs from the metadata list file
sampleList = {}
metadataFilename = '/projects/covidtrakr/tunc/biobot-SraRunTable.txt'
with open(metadataFilename) as infile:
    reader = csv.reader(infile, delimiter=",")
    next(reader)
    for row in reader:
        # site = row[41]
        state = row[18].split()[-1][0:4]
        sampleList[row[0]] = Sample(row[0], row[10].split('/')[0], state)

print('%d items found to process' % len(sampleList) )


# Analyse each report file 1 by 1
cumulatedDepth = []
cumulatedQuality = []
cumulatedIdx = []
for srr in sampleList.keys():
    tsvFilename = "/home/Tunc.Kayikcioglu/scratch/biobotOutput/%s/report/pos-coverage-quality.tsv" % srr
    depth = []
    quality = []
    with open(tsvFilename) as infile:
        reader = csv.reader(infile, delimiter="\t")
        for row in reader:
            depth.append(int(row[1]))
            quality.append(float(row[2]))
    
    #sampleList[srr].depth = depth            
    #sampleList[srr].quality = quality
    sampleList[srr].numPoorLoci = len( [ x for x in depth if x < 10 ] )
    sampleList[srr].numUncoveredLoci = len( [ x for x in depth if x < 1 ] )
    
    cumulatedDepth += depth
    cumulatedQuality += quality
    cumulatedIdx += list(range(1,len(quality)+1))
    
    if sampleList[srr].numPoorLoci < 1000:
        print('%s misses only %d loci.' % (srr,sampleList[srr].numPoorLoci) )



###############################################################################
plt.rcParams.update({'font.size': 14})
FDAblue = (0,124/255,186/255) # RGB color representation of the logo

SRRs = list(sampleList.keys())
sites = np.array([ sampleList[x].collectionSite for x in SRRs ])
numPoorLoci = np.array([ sampleList[x].numPoorLoci for x in SRRs ])
numBlindLoci = np.array([ sampleList[x].numUncoveredLoci for x in SRRs ])
dates = np.array([ sampleList[x].collectionDate for x in SRRs ])


# Generate a plot of poorly covered loci vs collection date
plt.plot(dates, numPoorLoci/1000, '.', color=FDAblue )
dateLims = plt.xlim()
dateTicks = np.arange(dateLims[0], dateLims[1], (dateLims[1]-dateLims[0])/10 )
plt.xticks(dateTicks, rotation=90)
plt.xlabel('Collection date')
plt.ylabel('# loci with < 10X coverage (1000)')

plt.savefig('./dateVSpoorLociCount.png', dpi=200, bbox_inches='tight')
plt.close()


# Generate a plot of uncovered loci vs collection date
plt.plot(dates, numBlindLoci/1000, '.', color=FDAblue )
plt.xticks(dateTicks, rotation=90)
plt.xlabel('Collection date')
plt.ylabel('# loci with no coverage (1000)')


plt.savefig('./dateVSblindLociCount.png', dpi=200,bbox_inches='tight')
plt.close()


# Generate a plot of poorly covered loci vs collection site
plt.plot(sites, numPoorLoci/1000, '.', color=FDAblue )
plt.xticks(rotation = 90, fontsize=8)
plt.xlabel('Collection site')
plt.ylabel('# loci with < 10X coverage (1000)')

plt.savefig('./siteVSpoorLociCount.png', dpi=200, bbox_inches='tight')
plt.close()


# Generate a plot of uncovered loci vs collection site
plt.plot(sites, numBlindLoci/1000, '.', color=FDAblue )
plt.xticks(rotation = 90, fontsize=8)
plt.xlabel('Collection site')
plt.ylabel('# loci with no coverage (1000)')

plt.savefig('./siteVSblindLociCount.png', dpi=200, bbox_inches='tight')
plt.close()



########################################################################
# 2D histograms for depth/quality vs genomic coordinate analyses
plt.rcParams['image.cmap'] = 'binary'


# Cumulated quality plot
plt.hist2d(cumulatedIdx, cumulatedQuality, bins=[np.arange(0,30000,50), np.arange(-1,45,1)], density=True)
plt.xlabel('Genome position (bp)')
plt.ylabel('Average read quality')
#plt.colorbar()

plt.savefig('./cumulatedQuality.png', dpi=200, bbox_inches='tight')
plt.close()



########################################################################
# Cumulated read depth plot
meanDepth = np.mean(cumulatedDepth)
depth2plot = [ min(x,2*meanDepth) for x in cumulatedDepth ]
plt.hist2d(cumulatedIdx, depth2plot, bins=[np.arange(0,30000,50), 50], density=True)
plt.xlabel('Genome position (bp)')
plt.ylabel('Average sequencing depth')
#plt.colorbar()

plt.savefig('./cumulatedDepth.png', dpi=200, bbox_inches='tight')
plt.close()



########################################################################
# 1D read depth histogram
FDAblue = (0,124/255,186/255) # RGB color representation of the logo

plt.hist(cumulatedDepth, bins=100, color=FDAblue, edgecolor=None, density=True )
plt.yticks([])
plt.xlabel('Average sequencing depth')
plt.ylabel('Relative frequency (cumulative)')

plt.savefig('./depthHistogram.png', dpi=200, bbox_inches='tight')
plt.close()



########################################################################
# 1D quality histogram
plt.hist(cumulatedQuality, bins=100, color=FDAblue, edgecolor=None, density=True )
plt.yticks([])
plt.xlabel('Average read quality')
plt.ylabel('Relative frequency (cumulative)')
plt.xlim(-1,45)

plt.savefig('./qualityHistogram.png', dpi=200, bbox_inches='tight')
plt.close()


print('%% with quality > 30: %.2f' % (100.0*len([ x for x in cumulatedQuality if x >= 30 ]) / len(cumulatedQuality)) )
print('%% with depth < 10X: %.2f' % (100.0*len([ x for x in cumulatedDepth if x < 10 ]) / len(cumulatedDepth)) )
