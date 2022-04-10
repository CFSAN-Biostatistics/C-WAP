import sys
import vcf
import vcf.utils
import pandas
import os.path

# sites
output = sys.argv[1]
min_af = float(sys.argv[2])
vcf_files = sys.argv[3:]

# SAMPLE, CHROM, POS, REF, ALT, ADREF, ADALT, DP

muts = {}
for f in vcf_files:
	reader = vcf.Reader(open(f, 'rb'))
	for i, r in enumerate(reader):
		for j, mp in enumerate(r.ALT):
			for k, s in enumerate(reader.samples):
				ref, alt = vcf.utils.trim_common_suffix(str(r.REF), str(r.ALT[j]))
				key = (r.CHROM, r.POS, ref, alt)

				if key not in muts:
					muts[key] = {}
					muts[key]['rows'] = []
					muts[key]['max_af'] = 0
					
				muts[key]['rows'].append([s, r.CHROM, r.POS, ref, alt, r.samples[k]['AD'][0], r.samples[k]['AD'][j+1], r.samples[k]['DP']])
				if r.samples[k]['DP'] > 0:
					muts[key]['max_af'] = max(muts[key]['max_af'], r.samples[k]['AD'][j+1] / r.samples[k]['DP'])

out_rows = []
for k in sorted(muts.keys()):
	if muts[k]['max_af'] >= min_af:
		out_rows.extend(muts[k]['rows'])

df = pandas.DataFrame(out_rows, columns=['sample','chrom','pos','ref','alt','adref','adalt','dp'])
df.to_csv(output, sep='\t', index=False)
