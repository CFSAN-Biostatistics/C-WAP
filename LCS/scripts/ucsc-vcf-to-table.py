import sys
import gzip
import re

h = gzip.open(sys.argv[1],'rt')
table = {}
for l in h:
	if l.startswith('#'):
		continue
	else:
		rows = l.rstrip().split("\t")
		pos = rows[1]
		ref = rows[3]
		alts = rows[4].split(",")
		for i,s in enumerate(rows[9:]):
			if int(s) == 0:
				mut = pos+":"+ref+">"+ref
			else:
				mut = pos+":"+ref+">"+alts[int(s)-1]
			if pos not in table:
				table[pos] = {}
			if mut not in table[pos]:
				table[pos][mut] = 0
			table[pos][mut] += 1

print("\t".join(["chrom","pos","ref","alt","adref","adalt","dp"]))
for pos in table.keys():
	dp = 0
	ads = {}
	ref = ""
	for mut in table[pos].keys():
		dp += table[pos][mut]
		change = re.sub(".*:","",mut)
		ref = re.sub(">.*","",change)
		alt = re.sub(".*>","",change)
		
		if alt not in ads:
			ads[alt] = 0
		ads[alt] += table[pos][mut]

	adref = 0
	if ref in ads:
		adref = ads[ref]
	for i in ads.keys():
		if i == ref:
			continue
		adalt = ads[i]
		print("\t".join(["NC_045512.2", pos, ref, i, str(adref), str(adalt), str(dp)]))

