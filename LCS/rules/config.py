
REF='data/refs/NC_045512.2.fasta'
VARIANT_GROUPS="data/variant_groups.tsv"
SITES_MIN_AF=0.2

# gisaid data: https://www.gisaid.org/
GISAID_FASTA='data/gisaid.fa.gz'

# pango designations
PANGO_DESIGNATIONS_VERSION='v1.2.60'
PANGO_DESIGNATIONS_REPO='https://github.com/cov-lineages/pango-designation'
PANGO_DESIGNATIONS_DIR='data/pango-designation'
PANGO_DESIGNATIONS='data/pango-designation/lineages.csv'

# covTree data
PB_URL='https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/'
PB_VERSION='2021-08-19'
PB='data/ucsc-sars-cov-2/public-{}.all.masked.pb.gz'.format(PB_VERSION)
PB_METADATA='data/ucsc-sars-cov-2/public-{}.metadata.tsv.gz'.format(PB_VERSION)

# primer list for amplicon trimming
PRIMERS_FA='data/artic_v3_primers.fa'
primer_trimming = 'primer_trimming' in config and config['primer_trimming']

import os.path

if 'markers' not in config:
	raise Exception("set you must set --config markers={pango|ucsc} argument")

if 'dataset' not in config:
	raise Exception("set you must set --config dataset=<mydataset> argument")

def get_tags_file(x):
	file = 'data/tags_pool_{}'.format(x)
	if not os.path.exists(file):
		raise Exception("tags file not found: "+file)
	return file

def get_tags(x):
	return open(get_tags_file(x)).read().rstrip().split("\n")

def get_variant_groups():
	return list(map(lambda l: l.split("\t")[0], open(VARIANT_GROUPS).read().rstrip().split("\n")[1:]))

def get_lineages(c):
	return list(map(lambda l: l.split("\t")[1], filter(lambda l: l.startswith(c+"\t"), open(VARIANT_GROUPS).read().rstrip().split("\n")[1:])))[0]

def get_target_sites():
	return 'outputs/lineage_sites/{}_sites.vcf.gz'.format(config['markers'])

def get_target_markers_table():
	return 'outputs/variants_table/{}-markers-table.tsv'.format(config['markers'])

def get_target_pools_table():
	return 'outputs/variants_table/pool_samples_{}.tsv'.format(config['dataset'])


