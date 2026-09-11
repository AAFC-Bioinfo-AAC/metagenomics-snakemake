'''
    Filename: amr_short_reads.smk
    Author: Katherine James-Gzyl and Devin Holman
    Date created: 2026/09/11
    Snakemake version: 9.20.0
    Python version: 3.8
'''
# -------------------------------------------------------------------
# CARD database setup
# -------------------------------------------------------------------

def project_abspath(path):
    """Resolve a path relative to PROJECT_ROOT unless it is already absolute."""
    path = os.path.expandvars(os.path.expanduser(os.fspath(path)))
    if os.path.isabs(path):
        return os.path.normpath(path)
    return os.path.abspath(os.path.join(os.fspath(PROJECT_ROOT), path))


def resolve_card_db():
    """Return the configured, preloaded CARD/RGI local database directory."""
    env_card = os.getenv("RGI_CARD", "").strip()
    cfg_card = config.get("card_latest")

    if env_card:
        card_db = env_card
    elif isinstance(cfg_card, str) and cfg_card.strip():
        card_db = cfg_card.strip()
    else:
        raise ValueError(
            "Set RGI_CARD in .env or card_latest in config.yaml before "
            "requesting an AMR target."
        )

    return project_abspath(card_db)


def card_db_file(relative_path):
    """Return a lazy Snakemake input function for a CARD database file."""
    def _get_card_db_file(wildcards):
        return os.path.join(resolve_card_db(), relative_path)

    return _get_card_db_file


# -------------------------------------------------------------------
# Validate the preloaded CARD database and KMA index
# -------------------------------------------------------------------

rule rgi_validate_database:
    input:
        card_json = card_db_file("card.json"),
        card_fasta = card_db_file("card_reference.fasta"),
        loaded_databases = card_db_file("loaded_databases.json"),
        kma_comp = card_db_file("bwt/card_reference/kma.comp.b"),
        kma_length = card_db_file("bwt/card_reference/kma.length.b"),
        kma_name = card_db_file("bwt/card_reference/kma.name"),
        kma_seq = card_db_file("bwt/card_reference/kma.seq.b")
    output:
        marker = f"{LOG_DIR}/rgi_card_db.validated"
    log:
        f"{LOG_DIR}/rgi/rgi_validate_database.log"
    conda:
        "../envs/rgi.yaml"
    params:
        card_db = lambda wildcards: resolve_card_db()
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})"
        mkdir -p "$(dirname {output.marker:q})"
        : > {log:q}

        card_db={params.card_db:q}

        if [[ ! -d "$card_db" ]]; then
            echo "ERROR: CARD database directory not found: $card_db" >> {log:q}
            exit 1
        fi

        for required_file in {input:q}; do
            if [[ ! -s "$required_file" ]]; then
                echo "ERROR: Required CARD/RGI database file is missing or empty: $required_file" >> {log:q}
                exit 1
            fi
        done

        tmpbase="${{TMPDIR:-/tmp}}"
        validation_dir="$(mktemp -d "$tmpbase/rgi_validate_XXXXXX")" || {{
            echo "ERROR: Could not create temporary validation directory" >> {log:q}
            exit 1
        }}

        cleanup() {{
            if [[ -n "${{validation_dir:-}}" && -d "$validation_dir" ]]; then
                rm -rf -- "$validation_dir"
            fi
        }}
        trap cleanup EXIT

        # RGI --local expects a directory named localDB beneath its working
        # directory. Use a temporary symlink so the project root is not modified.
        ln -s "$card_db" "$validation_dir/localDB"

        (
            cd "$validation_dir"
            rgi database --version --local
        ) >> {log:q} 2>&1

        echo "Validated CARD/RGI database: $card_db" >> {log:q}
        touch {output.marker:q}
        """


# -------------------------------------------------------------------
# Map cleaned paired-end reads to CARD with RGI/KMA
# -------------------------------------------------------------------

rule rgi_bwt:
    wildcard_constraints:
        sample = "[^/]+"
    input:
        R1 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz",
        R2 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz",
        database_marker = f"{LOG_DIR}/rgi_card_db.validated"
    output:
        # Intermediate files retained only when required by another rule.
        json = temp(
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/{{sample}}_paired.allele_mapping_data.json"
        ),
        bam = temp(
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/{{sample}}_paired.sorted.length_100.bam"
        ),
        bai = temp(
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/{{sample}}_paired.sorted.length_100.bam.bai"
        ),

        allele = (
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/"
            f"{{sample}}_paired.allele_mapping_data.txt"
        ),
        gene = (
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/"
            f"{{sample}}_paired.gene_mapping_data.txt"
        ),
        artifacts_stats = (
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/"
            f"{{sample}}_paired.artifacts_mapping_stats.txt"
        ),
        overall_stats = (
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/"
            f"{{sample}}_paired.overall_mapping_stats.txt"
        ),
        reference_stats = (
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/"
            f"{{sample}}_paired.reference_mapping_stats.txt"
        )
    params:
        outprefix = lambda wildcards: project_abspath(
            f"{CARD_RGI_OUTPUT_DIR}/{wildcards.sample}/{wildcards.sample}_paired"
        ),
        card_db = lambda wildcards: resolve_card_db()
    log:
        f"{LOG_DIR}/rgi/bwt_{{sample}}.log"
    threads:
        config.get("rgi_bwt", {}).get("threads", 4)
    conda:
        "../envs/rgi.yaml"
    shell:
        r"""
        set -euo pipefail

        # Resolve read and log paths before changing into temporary storage.
        workflow_dir="$PWD"
        read1={input.R1:q}
        read2={input.R2:q}
        log_file={log:q}

        [[ "$read1" = /* ]] || read1="$workflow_dir/$read1"
        [[ "$read2" = /* ]] || read2="$workflow_dir/$read2"
        [[ "$log_file" = /* ]] || log_file="$workflow_dir/$log_file"

        card_db={params.card_db:q}
        outprefix={params.outprefix:q}

        mkdir -p "$(dirname "$log_file")"
        mkdir -p "$(dirname "$outprefix")"
        : > "$log_file"

        tmpbase="${{TMPDIR:-/tmp}}"
        run_dir="$(mktemp -d "$tmpbase/rgi_bwt_XXXXXX")" || {{
            echo "ERROR: Could not create temporary RGI working directory" >> "$log_file"
            exit 1
        }}

        cleanup() {{
            if [[ -n "${{run_dir:-}}" && -d "$run_dir" ]]; then
                rm -rf -- "$run_dir"
            fi
        }}
        trap cleanup EXIT

        ln -s "$card_db" "$run_dir/localDB"

        echo "CARD/RGI database: $card_db" >> "$log_file"
        echo "RGI working directory: $run_dir" >> "$log_file"

        (
            cd "$run_dir"
            rgi bwt \
                -1 "$read1" \
                -2 "$read2" \
                -a kma \
                -n {threads} \
                -o "$outprefix" \
                --local \
                --clean
        ) >> "$log_file" 2>&1

        # Snakemake checks for the declared paths. This additional test prevents
        # empty files from being accepted as successful RGI results.
        for result_file in {output:q}; do
            if [[ ! -s "$result_file" ]]; then
                echo "ERROR: Expected RGI output is missing or empty: $result_file" >> "$log_file"
                exit 1
            fi
        done
        """
