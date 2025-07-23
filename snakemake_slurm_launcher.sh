#!/bin/bash
#SBATCH --job-name=snk_metagenomics
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --cluster=gpsc8
#SBATCH --partition=standard
#SBATCH --account=aafc_aac
#SBATCH --mem=2000
#SBATCH --time=8:00:00

source /gpfs/fs7/aafc/common/miniforge/miniforge3/etc/profile.d/conda.sh

conda activate snakemake-9.6.0
export PATH="$PWD/bin:$PATH"

  snakemake \
    --profile /abs/path/to/project/code/profiles/slurm \
    --configfile /abs/path/to/project/code/config/config.yaml \
    --conda-prefix /abs/path/to/the/conda/env/ \
    --printshellcmds \
    --keep-going \
    --report /abs/path/to/project/report.html