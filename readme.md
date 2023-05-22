# C-WAP
# CFSAN Wastewater Analysis Pipeline

[Kayikcioglu T, Amirzadegan J, Rand H, Tesfaldet B, Timme RE, Pettengill JB. Performance of methods for SARS-CoV-2 variant detection and abundance estimation within mixed population samples. PeerJ. 2023 Jan 26;11:e14596. doi: 10.7717/peerj.14596. PMID: 36721781; PMCID: PMC9884472.](https://pubmed.ncbi.nlm.nih.gov/36721781/)



**Given the [project timeline](https://www.fda.gov/food/whole-genome-sequencing-wgs-program/wastewater-surveillance-sars-cov-2-variants), C-WAP will no longer be under active  development or maintenance come June 30, 2023. [Freyja](https://github.com/andersen-lab/Freyja) or [Kallisto](https://github.com/pachterlab/kallisto) (two of the tools C-WAP incorporates) may be of interest.  Thank you for joining us on our analytic journey.**



C-WAP is a Nextflow, Python, and bash-based bioinformatics pipeline for the analysis of either long-read (ONT or PacBio) or short-read (Illumina) whole genome sequencing
data of DNA extracted from wastewater. It was developed for SARS-CoV2 and its variants.

C-WAP was developed by the United States Food and Drug Administration, Center for Food Safety and Applied Nutrition.


### Introduction

The CFSAN Wastewater Analysis Pipeline uses a reference-based alignment to create a matrix of
SNPs for a given set of samples and estimate the perentage of SC2 variants in the sample 

The process includes the following:
1. Designating a reference and NGS data in fastq format
2. Alignment of reads to the reference via Bowtie2
3. Taxonomy check via kraken2
4. Processing of alignment results via samtools
5. Detection of variant positions with ivar
6. Determine composition of variants via kallisto, linear regression, kraken2/bracken and freyja
7. Generate an html and pdf formatted summary of results



### Dependencies

User provided:
* [Conda3](https://docs.conda.io/en/latest/miniconda.html)
* [NextFlow v21.12](https://github.com/nextflow-io/nextflow/releases/tag/v21.12.1-edge)

Auto-fetched by C-WAP:
* [ghostscript](https://www.ghostscript.com)
* [kraken2 v2.1.2 ](https://github.com/DerrickWood/kraken2)
* [bracken](https://github.com/jenniferlu717/Bracken)
* [samtools v1.13 ](https://github.com/samtools/)
* [bcftools](https://github.com/samtools/bcftools)
* [ivar](https://github.com/andersen-lab/ivar)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
* [minimap2 v2.22](https://github.com/lh3/minimap2)
* [entrez-direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
* [pangolin](https://github.com/cov-lineages/pangolin)
* [Freyja](https://github.com/andersen-lab/Freyja)
* [kallisto](https://github.com/pachterlab/kallisto)
* [wkhtmltopdf](https://github.com/wkhtmltopdf)

./startWorkflow.nf assumes that conda, nextflow and gs executables are available in the search path. All other dependencies are imported via conda. To trigger this installation, execute `c-wap/prepare_envs.sh` (no root privileges needed). The series of acquisition of the dependencies and creating of the env's might take a substantially long time (potentially hours). Analysis calls to C-WAP will make use of the cached env's stored under the c-wap/conda subdirectory and are expected to finish substantially faster. To update the dependencies (more often required for tools with a biological DB) delete the respective environement folder `c-wap/conda/env-<toolName>` and re-execute `c-wap/prepare_envs.sh`.


### Installation

Install nextflow and conda. Afterwards, download c-wap repository and save. In addition, also obtain a copy of the kraken2 standard DB (~50GB). You can either download a copy from here: https://benlangmead.github.io/aws-indexes/k2 or you can compile your own version following the instructions under the "Standard Kraken2 Database" section here: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown

Afterwards, set $K2_STD_DB_PATH to your download location.


### Usage 

The driver script is `startWorkflow.nf` and a standard execution with paired end illumina reads would be:  
`startWorkflow.nf --platform i --primers path/to/bed --in path/to/fastq/ --out path/to/outputDir`

Note that the input files need to be compressed fastq.gz formats only and .sam, .bam or uncompressed .fastq files result in a fatal runtime error.


### Output

C-WAP produces a number of files from the various processing steps.  

`sorted.stats` - Samtools stats output from aligned but untrimmed reads  
`kallisto.out` - Python-parsed summary of the kallisto lineage abundance estimates  
`deconvolution.output` - Output of linear deconvolution method for estimating variant composition  
`linearDeconvolution_abundance.csv` - Linear deconvolution estimates of variant composition  
`freyja.demix` - Lineage abundance estimate generated by Freyja  
`kallisto_abundance.tsv` - Kallisto estimates of variant composition  
`k2-majorCovid.out` - Covid-specific kraken2 output with major lineages identified, against majorCovid DB  
`k2-majorCovid_bracken.out` - Bracken lineage abundance estimates, against majorCovid DB  
`pangolin_lineage_report.csv` - Pangolin lineage prediction for the consensus sequence  
`consensus.fa` - consensus fasta file generated by ivar  
`calls.vcf.gz` - Variant call file generated by bcftools  
`pos-coverage-quality.tsv` - QC metrics on coverage and quality obtained from the pileup file  
`rawVarCalls.tsv` - Variant calls generated by iVar, vcg equivalent of samtools  
`k2-std.out` - Kraken2 output with standard database  
`absCounts.csv` - Absolute counts table of uncovered and undercovered genes  
`scaledCounts.csv` - Scaled counts table of uncovered and undercovered genes  
`report` - Standalone directory containing html and pdf summary report  


### Note about variant composition

Variant composition analyses should be interpreted with caution where they should be treated as suspect if there are substantial gaps in coverage across the reference genome and/or a lack of sequencing depth.  The linear deconvolution and kraken2/bracken covid method are internally developed methods and under testing and validation.  

### GISAID citation
We gratefully acknowledge all data contributors, i.e., the Authors and their Originating laboratories responsible for obtaining the specimens, and their Submitting laboratories for generating the genetic sequence and metadata and sharing via the GISAID Initiative, some of which this research utlizes.

Khare, S., et al (2021) GISAID’s Role in Pandemic Response. China CDC Weekly, 3(49): 1049-1051. doi: 10.46234/ccdcw2021.255 PMCID: 8668406

### Citing C-WAP
[Kayikcioglu T, Amirzadegan J, Rand H, Tesfaldet B, Timme RE, Pettengill JB. Performance of methods for SARS-CoV-2 variant detection and abundance estimation within mixed population samples. PeerJ. 2023 Jan 26;11:e14596. doi: 10.7717/peerj.14596. PMID: 36721781; PMCID: PMC9884472.](https://pubmed.ncbi.nlm.nih.gov/36721781/)

### License

See the LICENSE.txt file included in the C-WAP Pipeline distribution.

