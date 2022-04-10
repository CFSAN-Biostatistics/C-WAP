
include: 'rules/config.py'
include: 'rules/ucsc-matrix-pipe.py'
include: 'rules/pango-matrix-pipe.py'
include: 'rules/pools-pipe.py'

rule all:
	input: expand('outputs/decompose/{x}.out', x=config['dataset'])
