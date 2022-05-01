import sys
import vcf
import vcf.utils
import glob
import pandas

muts = set()
min_af = float(sys.argv[1])

df = pandas.read_csv(sys.argv[2], sep='\t')
df['af'] = df['adalt']/df['dp']
df = df[df.af>=min_af]

df = df[['chrom','pos','ref','alt']]
df = df.drop_duplicates()
df = df.sort_values(['chrom','pos','ref','alt'])

print('##fileformat=VCFv4.2')
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
for i,r in df.iterrows():
	print("\t".join([r[0],str(r[1]),'.',r[2],r[3],'.','.','.']))

