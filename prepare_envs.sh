#!/bin/bash -e

# Check if the conda directory exists. If not, make one.
if ! [[ -d ./conda ]]; then
    mkdir ./conda
fi

conda_flags="-c bioconda -c conda-forge -c defaults --mkdir --yes --quiet"

if ! [[ -d conda/env-bowtie2 ]]; then
    conda create ${conda_flags} --prefix conda/env-bowtie2 bowtie2
fi

if ! [[ -d conda/env-minimap2 ]]; then
    conda create ${conda_flags} --prefix conda/env-minimap2 minimap2
fi

if ! [[ -d conda/env-ivar ]]; then
    conda create ${conda_flags} --prefix conda/env-ivar ivar=1.3.1
fi

if ! [[ -d conda/env-samtools ]]; then
    conda create ${conda_flags} --prefix conda/env-samtools samtools
fi

if ! [[ -d conda/env-kraken2 ]]; then
    conda create ${conda_flags} --prefix conda/env-kraken2 kraken2=2.1.2 bracken=2.7
fi

if ! [[ -d conda/env-python ]]; then
    conda create ${conda_flags} --prefix conda/env-python matplotlib scikit-learn pandas
fi

if ! [[ -d conda/env-bcftools ]]; then
    conda create ${conda_flags} --prefix conda/env-bcftools bcftools
fi

if ! [[ -d conda/env-pangolin ]]; then
    conda create ${conda_flags} --prefix conda/env-pangolin pangolin=4.1.2
fi

if ! [[ -d conda/env-kallisto ]]; then
    conda create ${conda_flags} --prefix conda/env-kallisto kallisto
fi

if ! [[ -d conda/env-freyja ]]; then
    conda create ${conda_flags} --prefix conda/env-freyja freyja=1.3.8
fi

if ! [[ -d conda/env-entrez-direct ]]; then
    conda create ${conda_flags} --prefix conda/env-entrez-direct entrez-direct
fi

if ! [[ -d conda/env-gs-wkhtmltopdf ]]; then
    conda create ${conda_flags} --prefix conda/env-gs-wkhtmltopdf openssl=1.0 wkhtmltopdf=0.12.4 ghostscript=9.54
fi


exit 0

