import sys
import Bio.SeqIO
import gzip
from collections import OrderedDict
import os
import re

file_in = sys.argv[1]
vg_tsv = sys.argv[2]
pango_designation_csv = sys.argv[3]
out_prefix = sys.argv[4]

class LRU(OrderedDict):
    'Limit size, evicting the least recently looked-up key when full'

    def __init__(self, maxsize=128):
        self.maxsize = maxsize
        super().__init__()

    def __getitem__(self, key):
        value = super().__getitem__(key)
        self.move_to_end(key)
        return value

    def __setitem__(self, key, value):
        if key in self:
            self.move_to_end(key)
        super().__setitem__(key, value)
        if len(self) > self.maxsize:
            oldest = next(iter(self))
            del self[oldest]

def lineage_file(lineage):
    return out_prefix+"/"+lineage+".fa.gz"

def open_lineage_file(lineage):
    return gzip.open(lineage_file(lineage), 'at')

fh_cache = LRU(maxsize=500) # 500 file handles

vgs = {}
h = open(vg_tsv, 'rt')
for l in h.readlines():
	name,ls = l.rstrip().split("\t")
	ls = ls.split(",")
	for x in ls:
		vgs[x] = name
print(vgs)
lineages = {}
h = open(pango_designation_csv,"rt")
for l in h:
    fn,lineage = l.rstrip().split(",")
    if fn == 'taxon':
        continue # header
    lineages[fn] = lineage

h = gzip.open(file_in,"rt")
for i, r in enumerate(Bio.SeqIO.parse(h, "fasta")):
	m = re.sub("^hCoV-19/","",re.sub("\\|.*","",r.id))
	if m in lineages:
		lin = lineages[m]
		if lin in vgs:
			cons = vgs[lin]
			if cons not in fh_cache:
				fh_cache[cons] = open_lineage_file(cons)
			
			seq = str(r.seq.upper()) # upper-case
			seq = re.sub('[^ACGTN]', 'N', seq)  # remove ambigous bases
			fh_cache[cons].write((">"+r.id+"\n"+seq+"\n"))
