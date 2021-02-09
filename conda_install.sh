#!/bin/bash
set -ex

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Install tools, each in a separate  environment.
conda create -y -n parallel parallel
conda create -y -n quast quast==4.6.3
conda create -y -n busco busco=5.0.0
conda create -y -n bwa bwa==0.7.17
conda create -y -n samtools samtools==1.9.0
conda create -y -n mosdepth mosdepth==0.3.1
conda create -y -n bedtools bedtools==2.28.0
conda create -y -n ruby ruby==2.7.2 make gcc_linux-64
conda create -y -n r r-base==4.0.3 pandoc

# Install R libraries.
conda run -n r R -e "install.packages(c('ggplot2', \
    'ggpubr', 'rmarkdown', 'scales', 'tidyr'), \
    repos='https://cran.ma.imperial.ac.uk/')"

# Install Ruby libraries.
conda run -n ruby gem install numo-narray

# Clean up.
conda clean --all