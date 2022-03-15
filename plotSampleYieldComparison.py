#!/usr/bin/env python3

# Generates a summary plot comparing the relative number of covid reads in all samples
# Useful to judge if the SNR ove -ve control is sufficient.

import numpy as np
import matplotlib.pyplot as plt
import sys


dataIn = sys.argv[1:]
sample_names = dataIn[0::2]
numCovidReads = [float(x) for x in dataIn[1::2]]
numSamples = len(numCovidReads)


# Plot the absolute covid read counts as a bar graph
FDAblue = (0, 124/255, 186/255)  # RGB color representation of the logo
plt.rcParams.update({'font.size': 14})

plt.bar(np.arange(1,numSamples+1), numCovidReads, color=FDAblue)
plt.xlabel('Sample number')
plt.xticks(np.arange(1,numSamples+1, int(max(1,numSamples/10))))

# If read numbers are too high, scale the axes for a better view
ylocs, labels = plt.yticks()
if ylocs[-1] > 5000:
    plt.ylabel('# covid reads (1000)')
    plt.yticks(ylocs, (ylocs/1000).astype('int') )
else:
    plt.ylabel('# covid reads')


# SNR calculation
# Check if there is a negative control sample
def is_neg_control(name):
    for neg_word in ['water', 'negative', 'blank']:
        if neg_word in name.lower():
            return True
    return False
    

neg_control_idx = []
for name_idx in range(len(sample_names)):
    if is_neg_control (sample_names[name_idx]):
        neg_control_idx.append(name_idx)


# If it was clearly marked, then draw baseline levels as well.
if len(neg_control_idx)>0:
    mean_bg_reads = np.mean([ numCovidReads[x] for x in neg_control_idx ])
    plt.plot([0,numSamples+1], [mean_bg_reads,mean_bg_reads], 'k--')
    plt.xlim(0,numSamples+1)
    
    text_ylevel = max(1.2*mean_bg_reads, ylocs[-1]/3)
    for i in range(len(sample_names)):
        if i not in neg_control_idx:
            plt.text(i+1, text_ylevel, 'SNR=%.1f' % (numCovidReads[i]/mean_bg_reads), rotation=90,
                    horizontalalignment='center', verticalalignment='center')
        else:
            plt.text(i+1, text_ylevel, '-control', rotation=90, horizontalalignment='center',
                    verticalalignment='center')


plt.tight_layout()
plt.savefig('./covidReadsSummary.png', dpi=200)

