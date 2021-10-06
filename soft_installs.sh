#!/bin/bash
##
## Installs the soft stack for discerns_manuscript
##
## Izaskun Mallona
## GPLv3

# mkdir -p ~/virtenvs
# cd ~/virtenvs

# virtualenv -p python3 discerns

# source discerns/bin/activate
# which python
# which pip

# pip install numpy pandas
# pip install snakemake


# ## sra tools

# mkdir -p ~/soft/sra-tools
# cd $_

# ## precompiled for ubuntu i686
# wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.1/sratoolkit.2.11.1-ubuntu64.tar.gz
# tar xzvf sratoolkit.2.11.1-ubuntu64.tar.gz

# deactivate


## switch to conda (tmux condadiscerns):

source ~/miniconda3/bin/activate

conda create --name discerns_env

conda activate discerns_env

conda install -n discerns_env -c conda-forge mamba


mamba repoquery depends r-base -c conda-forge
mamba install r-base -c conda-forge
conda list | grep r-base

mamba repoquery depends star -c bioconda 
mamba install star -c bioconda

mamba repoquery depends rsem -c bioconda 
mamba install rsem -c bioconda
conda list | grep r-base ## still 4.1

mamba install bowtie2 -c bioconda
conda list | grep r-base ## went down to 4.0.5

mamba install hisat2 -c bioconda
mamba install sra-tools -c bioconda

# for tophat, python 2.7 is needed

conda create --name tophat_discerns_py27 python=2.7
conda activate tophat_discerns_py27

mamba install python=2.7 tophat -c bioconda -n tophat_discerns_py27
conda deactivate ## closing tophat_discerns_py27, so back to parent discerns_env

# EQP is still missing, I guess comes from https://github.com/Novartis/EQP-cluster
mamba repoquery depends EQP-QM -c bioconda
mamba install openjdk git -c conda-forge
mamba install bedtools -c bioconda

# then, clone EQP-QM and put it somewhere in the local path,
# https://github.com/novartis/EQP-QM (using git from the virtenv,
# plus java etc for exec)

# EQP end

conda env export > discerns_env.yaml
conda env export -n tophat_discerns_py27 > tophat_discerns_py27.yaml

mamba install -c conda-forge -c bioconda snakemake

## outer discerns_env
conda deactivate
