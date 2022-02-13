#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys


histValues = []
for line in sys.stdin:
    histValues.append(int(line.strip()))


plt.rcParams.update({'font.size': 14})
FDAblue = (0, 124/255, 186/255)  # RGB color representation of the logo

plt.hist(histValues, 100, color=FDAblue)
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

