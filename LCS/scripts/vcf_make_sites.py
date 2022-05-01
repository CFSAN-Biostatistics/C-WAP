import sys
import vcf
import vcf.utils
import glob
import pandas

muts = set()
min_af = float(sys.argv[1])

for f in sys.argv[2:]:
	reader = vcf.Reader(open(f, 'rb'))
	for i, r in enumerate(reader):
		for j, mp in enumerate(r.ALT):
			for k, s in enumerate(reader.samples):
				af = r.samples[k]['AD'][j+1] / r.samples[k]['DP']
				if af >= min_af:
					ref, alt = vcf.utils.trim_common_suffix(str(r.REF), str(r.ALT[j]))
					muts.add((r.CHROM, r.POS, ref, alt))

df = pandas.DataFrame(muts)
df = df.sort_values([0,1,2,3])
print('##fileformat=VCFv4.2')
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
for i,r in df.iterrows():
	print("\t".join([r[0],str(r[1]),'.',r[2],r[3],'.','.','.']))

