#! /usr/bin/env nextflow

// Input is grouped in subfolders to reduce the overhead due to thousands of slurm jobs taking 10s each.
inSubfolders = Channel.fromPath("$params.inpath/*", type: 'any', checkIfExists: true)


process freyja {
	time = '20 min'
	cpus = 2
	memory = '8 GB'
	errorStrategy = 'retry'
	maxRetries = 10
	
	input:
		path('subfolder') from inSubfolders
	
	output:
		path('demix/*') into outfiles
	
	publishDir "$params.outpath/", mode: 'copy', overwrite: true
	conda 'freyja=1.3.7 samtools=1.15'
	
	shell:
	"""
		unset PYTHONPATH
		mkdir ./demix
		for fileIdx in \$(ls subfolder/depths/* | awk -F '/' '{print \$NF}' | awk -F '.' '{print \$1}'); do
			freyja demix subfolder/variants/\$fileIdx.tsv subfolder/depths/\$fileIdx.tsv --output ./demix/\$fileIdx --confirmedonly
		done
	"""
}

