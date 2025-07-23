# Snakemake Profile for running on SLURM

## Directory layout 
```
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
- Must invoke with ` --profile /absolute/path/to/the/my_pipeline/profiles/slurm`
### Example Profile
`profiles/slurm/config.yaml`
```bash
###How Snakemake assign resources to rules ###
cores: 60 # total number of cores Snakemake can request at any time
jobs: 10 # max amount of jobs Snakemake runs at once
latency-wait: 60 # gives time for containers to start, and for other I/O delays
rerun-incomplete: true
quiet: false # Makes Snakemake output more verbose about its operations
retries: 2               # so jobs killed by IO hiccups are auto retried
max-jobs-per-second: 2  # slowing down submission rate (tune for your site!)
executor: slurm 

### Env Vars ###
envvars:
  TMPDIR: "/gpfs/fs7/aafc/scratch/$USER/tmpdir"

default-resources:
  - mem_mb=8000             
  - slurm_partition=standard  # SLURM partition to use by default
  - slurm_account=aafc_aac # SLURM account
  - slurm_cluster=gpsc8 #SLURM cluster
  - runtime=60       # minutes
  - mem_mb=4000 # default memory per job (change as desired)
  - cpus=1

### Env modules ###
# use-envmodules: false # use-envmodules: true only if: Your system disables Conda/containers and expects you to use module load bioinfo-tool for each step or you have a properly configured profile and know which modules are needed for every rule.

### Conda ###
use-conda: true
conda-frontend: mamba   # if you have mamba installed

### Resource scopes ###
# Affects how Snakemake defines resource limits
# Note that cores and threads are always considered local 
# set-resource-scopes:
    cores: local 

### Containerization if using containers like Docker ###
# use-singularity: true 
# singularity-prefix-dir: 
# singularity-args:
# cleanup-containers: false 

### Rule-specific configs ###
# These go into the cluster_config.yaml
## Per rule resources
set-resources:
  fastp_pe:
    cpus: 2
    mem_mb: 4000
    runtime: 40
    slurm_partition: standard
    slurm_account: aafc_aac

  bowtie2_align:
    cpus: 24
    mem_mb: 48000
    runtime: 30
    slurm_partition: standard
    slurm_account: aafc_aac

  extract_unmapped_fastq:
    cpus: 60
    mem_mb: 64000
    runtime: 30
    slurm_partition: standard
    slurm_account: aafc_aac

  sortmerna_pe:
    cpus: 48
    mem_mb: 32000
    runtime: 360
    slurm_partition: standard
    slurm_account: aafc_aac

  kraken2:
    cpus: 2
    mem_mb: 600000
    runtime: 30 
    slurm_partition: large
    slurm_account: aafc_aac__large

  bracken:
    cpus: 2
    mem_mb: 4000
    runtime: 10
    slurm_partition: standard
    slurm_account: aafc_aac

  combine_bracken_outputs:
    cpus: 1
    mem_mb: 2000
    runtime: 20
    slurm_partition: standard
    slurm_account: aafc_aac
   
  bracken_extract:
    cpus: 1
    mem_mb: 2000
    runtime: 10
    slurm_partition: standard
    slurm_account: aafc_aac

  rgi_reload_database:
    cpus: 1
    mem_mb: 2000
    runtime: 30
    slurm_partition: standard
    slurm_account: aafc_aac

  rgi_bwt:
    cpus: 20
    mem_mb: 64000
    runtime: 60
    slurm_partition: standard
    slurm_account: aafc_aac

  rna_spades:             
    cpus: 48
    mem_mb: 64000
    runtime: 240
    slurm_partition: standard
    slurm_account: aafc_aac

  rnaquast_busco:
    cpus: 4
    mem_mb: 16000
    runtime: 10
    slurm_partition: standard
    slurm_account: aafc_aac

  megahit_coassembly:
    cpus: 48
    mem_mb: 256000
    runtime: 240
    slurm_partition: standard
    slurm_account: aafc_aac

  index_coassembly:
    cpus: 8
    mem_mb: 16000
    runtime: 120
    slurm_partition: standard
    slurm_account: aafc_aac

  bowtie2_map_transcripts:
    cpus: 16
    mem_mb: 32000
    runtime: 30
    slurm_partition: standard
    slurm_account: aafc_aac

  assembly_stats_depth:
    cpus: 2
    mem_mb: 2000
    runtime: 30
    slurm_partition: standard
    slurm_account: aafc_aac

  prodigal_genes:
    cpus: 1
    mem_mb: 2000
    runtime: 60
    slurm_partition: standard
    slurm_account: aafc_aac

  featurecounts:
    cpus: 4
    mem_mb: 8000
    runtime: 10
    slurm_partition: standard
    slurm_account: aafc_aac