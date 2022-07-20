#!/usr/bin/env python3

# ONT only
# Plots the DNA flux (i.e. read count) vs. acquisition time
# with the help of timestamps reported by the instrument in 
# the fastq header line of each read.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import sys


if len(sys.argv) != 3:
    raise Exception('Incorrect call to the script.')
    
# Data regarding the current sample passed by UNIX
timestamp_filename = sys.argv[1]
outfilename = sys.argv[2]


df = pd.read_csv(timestamp_filename, header=None)
timestamps = [mdates.date2num(np.datetime64(x)) for x in list(df[0])]

locator = mdates.AutoDateLocator()
formatter = mdates.ConciseDateFormatter(locator)

plt.rcParams.update({'font.size': 14})
FDAblue = (0,124/255,186/255) # RGB color representation of the logo
plt.rcParams['image.cmap'] = 'binary'


# Generate a plot of the number of the reads mapping to reference vs time
plt.hist(timestamps, bins=50, color=FDAblue )
dateLims = plt.xlim()
dateTicks = np.arange(dateLims[0], dateLims[1], (dateLims[1]-dateLims[0])/10 )
plt.xticks(dateTicks, rotation=90)
plt.xlabel('Passage start time')
plt.ylabel('# Reads (1000)')
locs, labels = plt.yticks()
plt.yticks(locs, [int(int(x)/1000) for x in locs ])

plt.gca().xaxis.set_major_locator(locator)
plt.gca().xaxis.set_major_formatter(formatter)
plt.savefig(outfilename, dpi=200, bbox_inches='tight')
plt.close()
    
    