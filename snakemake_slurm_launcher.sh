#!/bin/bash
#SBATCH --job-name=mag_test_run
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --cluster=gpsc8
#SBATCH --partition=standard
#SBATCH --account=aafc_aac
#SBATCH --mem=2000
#SBATCH --time=288:00:00

source /gpfs/fs7/aafc/common/miniforge/miniforge3/etc/profile.d/conda.sh

conda activate snakemake-9.9.0
export PATH="$PWD/bin:$PATH"

  snakemake \
    --profile /gpfs/fs7/aafc/projects/J-003518_afc_rdar_malate/code/metagenomics-snakemake/profiles/slurm \
    --configfile /gpfs/fs7/aafc/projects/J-003518_afc_rdar_malate/code/metagenomics-snakemake/config/config.yaml \
    --conda-prefix /gpfs/fs7/aafc/projects/J-003518_afc_rdar_malate/code/malate-smk-conda-env \
    --printshellcmds \
    --latency-wait 120 