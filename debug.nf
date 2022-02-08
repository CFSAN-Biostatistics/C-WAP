#! /usr/bin/env nextflow

X = Channel.from(1,2,3,4)

process a {
	input:
		val(x) from X
	
	output:
		tuple val(y), val(x) into Y1, Y2
		
	shell:
	"""
		echo \$x
		y=\$x
	"""
}


process b {
	input:
		tuple val(y), val(x) from Y1
	
	output:
		tuple val(y), val(x) into Z
		
	shell:
	"""
		echo \$x
		y=\$x
	"""
}


process c {
	input:
		tuple val(y), val(x) from Z
		tuple val(y), val(x) from Y2
	
	output:
		val(y)
		
	shell:
	"""
		echo \$x
	"""
}


