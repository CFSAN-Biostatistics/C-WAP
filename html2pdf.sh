#!/bin/bash


# Generates a pdf version of the report file by rendering the html file
# The program was obtained from: https://wkhtmltopdf.org/downloads.html
# And unpacked locally by:
# rpm2cpio wkhtmltox-0.12.6-1.centos7.x86_64.rpm  | cpio -i -d

outDir=$1
sampleName=$2

echo Generating pdf output...
./wkhtmltopdf --enable-local-file-access --page-size Letter --margin-top 10mm --margin-bottom 0 \
	--margin-left 0 --margin-right 0 --print-media-type --title "C-WAP report" \
	$outDir/${sampleName}_report/report.html $outDir/${sampleName}_report/$sampleName.pdf
	
	