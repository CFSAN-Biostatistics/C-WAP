#! /bin/env python3

# A custom script to obtain a fastq file from a sam file via CIGAR strings
# This is useful when there is a bam file with soft clipping
# that needs to be converted into a fastq file

import csv, sys, re


if len(sys.argv) != 3:
    raise Exception('Incorrect call to the script.')
    
samfilename = sys.argv[1]
fastqfilename = sys.argv[2]


with open(samfilename, 'r') as infile, open(fastqfilename, 'w') as outfile:
    reader = csv.reader(infile, delimiter="\t", quoting=csv.QUOTE_NONE)
    read_counter=1
    for row in reader:
        if row[0][0] != '@':
            CIGAR = row[5]
            seq = row[9]
            quality = row[10]
            
            # Process the CIGAR string
            CIGAR_letters = [x for x in re.split('\d', CIGAR) if x]
            CIGAR_numbers = [x for x in re.split('\D', CIGAR) if x]
            
            if CIGAR_letters[0] == 'S':
                start_pos = int(CIGAR_numbers[0])
            else:
                start_pos = 0
            
            if CIGAR_letters[-1] == 'S':
                end_pos = -int(CIGAR_numbers[-1])
            else:
                end_pos = len(seq)
    
            trimmed_seq = seq[start_pos:end_pos]
            trimmed_quality = quality[start_pos:end_pos]
            
            # Output to the fastq file
            outfile.writelines(["@read%d\n" % read_counter, trimmed_seq, "\n+\n", trimmed_quality, "\n"])
            read_counter += 1


