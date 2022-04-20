#! /usr/bin/env python3
# Parses the bed file to detect uncovered genomic coordinates

import csv

def findUncoveredCoordinates (bedfile, is_trimmed):
    # Import the bed file and find the
    fwd_primers_start_pos = []
    rev_primers_end_pos = []

    with open(bedfile) as infile:
        reader = csv.reader(infile, delimiter="\t", quoting=csv.QUOTE_NONE)
        for row in reader:
            # Skip the empty rows in the bed file, if any
            if len(row) < 6:
                continue
            
            start = int(row[1])
            end = int(row[2])
            strand = row[5]
            
            # If the primers are going to be trimmed out, then the first and last primers' 
            # footprints are also uncovered.
            if is_trimmed:
                trim_length = end - start
            else:
                trim_length = 0

            if strand == '+':
                fwd_primers_start_pos.append(start + trim_length)
            else:
                rev_primers_end_pos.append(end-1-trim_length)


    # Uncovered region is defined by the genome compartment before the first fwd primer
    # or beyond the rev primer
    gaps = []
    if len(fwd_primers_start_pos) > 0:
        gaps.append((1,min(fwd_primers_start_pos)))
        
    if len(rev_primers_end_pos) > 0:
        gaps.append((max(rev_primers_end_pos),29903))
        
    return gaps



if __name__ == '__main__':
    import sys
    bedfile = sys.argv[1]
    is_trimmed = bool(int(sys.argv[2]))
    
    # Print the uncovered position's coordinate bounds to stdout
    limits = findUncoveredCoordinates(bedfile, is_trimmed)
    for (a,b) in limits:
        print('%d %d' % (a,b))
    
    