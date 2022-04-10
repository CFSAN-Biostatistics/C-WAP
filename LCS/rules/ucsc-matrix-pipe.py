
import gzip
import pandas

rule download_pb:
	output: PB
	params: full_url=os.path.join(PB_URL,(PB_VERSION.replace('-','/')), os.path.basename(PB))
	shell: 'wget {params.full_url} -O {output}'

rule download_metadata:
	output: PB_METADATA
	params: full_url=os.path.join(PB_URL,(PB_VERSION.replace('-','/')), os.path.basename(PB_METADATA))
	shell: 'wget {params.full_url} -O {output}'

rule sample_list:
	input: PB_METADATA
	output: temp('outputs/ucsc-vcf/{x}.sample-list')
	params: lin=lambda wc: get_lineages(wc.x)
	run:
		h = pandas.read_csv(input[0], sep='\t')
		out = open(output[0], 'wt')
		ls = params.lin.split(",")
		rows = h.loc[h['pangolin_lineage'].isin(ls)]
		out.write("\n".join(rows['strain'])+"\n")

rule lineage_vcf:
	input: pb=PB, l='outputs/ucsc-vcf/{x}.sample-list'
	output: vcf='outputs/ucsc-vcf/{x}.vcf.gz'
	shell: 'matUtils extract -i {input.pb} -s {input.l} -a 3 -v {output}'

rule ucsc_variants_to_table:
	input: 'outputs/ucsc-vcf/{x}.vcf.gz'
	output: 'outputs/variants_table/ucsc/{x}.tsv'
	shell: 'python scripts/ucsc-vcf-to-table.py {input} > {output}'

rule ucsc_gather_tables:
	input: expand('outputs/variants_table/ucsc/{x}.tsv', x=get_variant_groups())
	output: 'outputs/variants_table/ucsc-markers-table.tsv'
	run:
		ds = []
		for i,s in enumerate(get_variant_groups()):
			d = pandas.read_csv(input[i], sep='\t')
			d['sample'] = s
			ds.append(d)
		ds = pandas.concat(ds)
		ds.to_csv(output[0], sep='\t', index=False)

rule ucsc_make_sites:
	input: 'outputs/variants_table/ucsc-markers-table.tsv'
	output: 'outputs/lineage_sites/ucsc_sites.vcf.gz'
	shell: 'python scripts/ucsc_make_sites.py {SITES_MIN_AF} {input} | bgzip -c > {output} && tabix {output}'
