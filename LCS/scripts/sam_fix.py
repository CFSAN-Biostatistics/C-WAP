
import sys
import gzip
import pysam

# add fake base qualities 'IIII'
# add RG tag

sam_input = sys.argv[1]
rgtag = sys.argv[2]
sam_output = sys.argv[3]

sf = pysam.AlignmentFile(sam_input, "rb")

header = sf.header.as_dict()
header['RG'] = [{ 'ID': rgtag, 'SM': rgtag, 'LB': rgtag }]

sf_out = pysam.AlignmentFile(sam_output, "wb", header=header)
for r in sf:
	# add qualities
	r.query_qualities = [40] * len(r.query_sequence)
	r.set_tag('RG',rgtag)
	sf_out.write(r)

