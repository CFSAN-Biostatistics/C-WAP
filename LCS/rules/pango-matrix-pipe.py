
rule repo:
	output: PANGO_DESIGNATIONS, directory(PANGO_DESIGNATIONS_DIR)
	shell: 'git clone --depth 1 --branch {PANGO_DESIGNATIONS_VERSION} {PANGO_DESIGNATIONS_REPO} {PANGO_DESIGNATIONS_DIR}'

rule split_fasta_lineages:
	input: fa=GISAID_FASTA, desig=PANGO_DESIGNATIONS
	output: expand('outputs/fasta_lineages/{x}.fa.gz', x=get_variant_groups())
	params: prefix='outputs/fasta_lineages/'
	shell: 'python scripts/split_fasta_by_lineages2.py {input.fa} {VARIANT_GROUPS} {input.desig} {params.prefix}'

rule minimap2:
	input: 'outputs/fasta_lineages/{x}.fa.gz'
	output: temp('outputs/map/raw/{x}.bam')
	threads: 8
	resources: mem_gb=8
	shell: 'minimap2 -t {threads} -ax map-pb --for-only --secondary=no {REF} {input} | samtools view -F256 -F4 -F2048 -h | picard -Xmx{resources.mem_gb}G SortSam I=/dev/stdin O={output} SO=coordinate'

rule sam_fix:
	input: 'outputs/map/raw/{x}.bam'
	output: 'outputs/map/fix/{x}.bam'
	params: rg='{x}'
	shell: 'python scripts/sam_fix.py {input} {params.rg} {output} && samtools index {output}'

rule pango_mutect:
	input: 'outputs/map/fix/{x}.bam'
	output: 'outputs/mutect/{x}.vcf.gz'
	resources: mem_gb=8
	log: 'outputs/mutect/{x}.log'
	params: jobname='pango_mutect_p1'
	shell: '''gatk --java-options '-Xmx{resources.mem_gb}G' Mutect2 -R {REF} -O {output} -I {input} --max-reads-per-alignment-start 200 -mnp-dist 0 >{log} 2>&1'''

rule pango_make_sites1:
	input: expand('outputs/mutect/{x}.vcf.gz', x=get_variant_groups())
	output: temp('outputs/lineage_sites/pango_sites.temp.vcf.gz')
	params: jobname='pango_make_sites1'
	resources: mem_gb=1
	shell: 'python scripts/vcf_make_sites.py {SITES_MIN_AF} {input} | bgzip -c > {output} && tabix {output}'

rule pango_mutect_part2:
	input: bam='outputs/map/fix/{x}.bam', sites='outputs/lineage_sites/pango_sites.temp.vcf.gz'
	output: 'outputs/mutect_v2/{x}.vcf.gz'
	resources: mem_gb=8
	params: jobname='pango_mutect_p2'
	log: 'outputs/mutect_v2/{x}.log'
	shell: '''gatk --java-options '-Xmx{resources.mem_gb}G' Mutect2 -R {REF} -O {output} -I {input.bam} --alleles {input.sites} --max-reads-per-alignment-start 200 -mnp-dist 0 >{log} 2>&1'''

rule pango_vcf_to_table:
	input: vcfs=expand('outputs/mutect_v2/{x}.vcf.gz',x=get_variant_groups())
	output: tbl='outputs/variants_table/pango-markers-table.tsv'
	params: jobname='pango_vcf_to_table'
	resources: mem_gb=1
	shell: 'python scripts/pango-variants_to_table.py {output.tbl} {SITES_MIN_AF} {input.vcfs}'

rule pango_make_sites2:
	input: 'outputs/variants_table/pango-markers-table.tsv'
	output: 'outputs/lineage_sites/pango_sites.vcf.gz'
	params: jobname='pango_make_sites2'
	resources: mem_gb=1
	shell: 'python scripts/ucsc_make_sites.py {SITES_MIN_AF} {input} | bgzip -c > {output} && tabix {output}'

