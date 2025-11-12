# Snakemake Profile for running on SLURM

## Directory layout

``` bash
metagenomics_pipeline/
├── Workflow/
│   └── Snakemake
    └── logs
    └── scripts
    └── envs
      └── tool1.yaml
      └── too2.yaml
      └── tool3.yaml
      └── ... 
├──resources 
│ 
├── profiles/
│   └── slurm/
│       └── config.yaml         ← profile config
├── config/
│   └── config.yaml             ← workflow data/sample config
├── run_snakemake.sh            ← your SLURM launcher
└── ...                         
```

## Main profile

- Controls how Snakemake schedules jobs on SLURM
- Must invoke with `--profile /absolute/path/to/the/my_pipeline/profiles/slurm`

### Example Profile

`profiles/slurm/config.yaml`

```bash
## For Snakemake to recognize this file it must be named config.yaml.
## Specific information to our HPC configuration has been replaced in this example with <CAPITAL LETTERS>
## Please edit to include the configuration of the cluster you are using

### How Snakemake will run on SBATCH ###
cores: 60
jobs: 10 
latency-wait: 60 
rerun-incomplete: true
retries: 2              
max-jobs-per-second: 2 
executor: slurm

# Prevent rerunning jobs just for Snakefile edits
## flags available [input, mtime, params, software-env, code, resources, none]
rerun-triggers: [input, params, software-env]

### Env Vars ###
envvars:
  TMPDIR: "/<PATH>/<TO>/<SCRATCH>/${USER}/tmpdir"

default-resources:
  - slurm_account=<ACCOUNT_NAME>
  - slurm_partition=<PARTITION_NAME>
  - slurm_cluster=<CLUSTER_NAME>
  - slurm_qos=<QOS_LEVEL>      # e.g., 'low' if jobs are held in queue for long
  - runtime=<RUNTIME_MINUTES>  # e.g., 60
  - mem_mb=<MEMORY_MB>         # e.g., 4000

### Env modules ###
# use-envmodules: false 

### Conda ###
use-conda: true
conda-frontend: mamba   

### Resource scopes ###
set-resource-scopes:
  cores: local 

# Reusable Slurm Blocks (anchors)
# Standard partition/account/cluster used by most rules
_slurm_std: &slurm_std
  slurm_partition: <PARTITION_NAME>
  slurm_account: <ACCOUNT_NAME_standard> # e.g., standard, large memory 
  slurm_cluster: <CLUSTER_NAME>

# Large memory partition/account/cluster used by some rules
_slurm_large: &slurm_large
  slurm_partition: <PARTITION_NAME>
  slurm_account: <ACCOUNT_NAME_large> # e.g., standard, large memory 
  slurm_cluster: <CLUSTER_NAME>

## Per rule resources
set-resources:
  fastp_pe:
    <<: *slurm_std
    mem_mb: 4000
    runtime: 40

  bowtie2_align:
    <<: *slurm_std
    mem_mb: 24000
    runtime: 90

  extract_unmapped_fastq:
    <<: *slurm_std
    mem_mb: 18000
    runtime: 60

  kraken2:
    <<: *slurm_large
    mem_mb: 840000
    runtime: 120
