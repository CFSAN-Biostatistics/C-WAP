#!/usr/bin/env python3

# Generates a summary plot comparing the relative number of covid reads in all samples
# Useful to judge if the SNR ove -ve control is sufficient.

import numpy as np
import matplotlib.pyplot as plt
import sys


dataIn = np.array(sys.argv[1:], dtype='float')
numTotalReads = dataIn[0::2]
pctCovidReads = dataIn[1::2]
numSamples = len(numTotalReads)
covidReads = numTotalReads*pctCovidReads/100


FDAblue = (0, 124/255, 186/255)  # RGB color representation of the logo
plt.rcParams.update({'font.size': 14})

plt.bar(np.arange(1,numSamples+1), covidReads, color=FDAblue)
plt.xlabel('Sample number')
plt.xticks(np.arange(1,numSamples+1), np.arange(1,numSamples+1))

# If read numbers are too high, scale the axes for a better view
locs, labels = plt.yticks()
if locs[-1] > 2000:
    plt.ylabel('# covid reads (1000)')
    plt.yticks(locs, (locs/1000).astype('int') )
else:
    plt.ylabel('# covid reads')

plt.tight_layout()
plt.savefig('./covidReadsSummary.png', dpi=200)

