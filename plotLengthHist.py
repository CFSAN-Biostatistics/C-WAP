# Plots a histogram of all reads in a provided list of integers (i.e. lengths)

import numpy as np
import matplotlib.pyplot as plt
import sys


histValues = []
for line in sys.stdin:
    histValues.append(int(line.strip()))


plt.rcParams.update({'font.size': 14})
FDAblue = (0, 124/255, 186/255)  # RGB color representation of the logo

if max(histValues) > 500:
    plt.hist(histValues, np.arange(0,2000,20), color=FDAblue)
    plt.xlim(0,2000)
else:
    plt.hist(histValues, np.arange(0,1000,20), color=FDAblue)
    plt.xlim(0,1000)
plt.xlabel('Read length')


locs, labels = plt.yticks()
if locs[-1] > 2000:
    plt.ylabel('Number of reads (1000)')
    plt.yticks(locs, (locs/1000).astype('int') )
else:
    plt.ylabel('Number of reads')

plt.tight_layout()
plt.savefig('./readLengthHist.png', dpi=200)
plt.close()

