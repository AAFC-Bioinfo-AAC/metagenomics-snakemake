#!/bin/bash
#SBATCH --job-name=run_snakemake.sh
#SBATCH --output=run_snakemake_%j.out 
#SBATCH --error=run_snakemake_%j.err 
#SBATCH --cluster=<CLUSTER_NAME>
#SBATCH --partition=<PARTITION_NAME>
#SBATCH --account=<ACCOUNT_NAME>
#SBATCH --mem=<MEMORY_MB>         # e.g., 2000
#SBATCH --time=<HH:MM:SS>         # Must be long enough for completion of workflow 

source path/to/source/conda/common/miniforge/miniforge3/etc/profile.d/conda.sh

conda activate snakemake_env
export PATH="$PWD/bin:$PATH"

  snakemake \
    --profile absolute/path/to/profiles/slurm \
    --configfile absolute/path/to/config/config.yaml \
    --conda-prefix absolute/path/to/common/conda/metatranscriptomics-snakemake-conda \
    --printshellcmds \
    --latency-wait 120 