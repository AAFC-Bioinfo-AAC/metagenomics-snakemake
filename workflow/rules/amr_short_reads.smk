'''
    Filename: amr_short_reads.smk
    Author: Katherine James-Gzyl
    Date created: 2025/07/16
    Snakemake version: 9.9.0
    Python version: 3.8
'''

CARD_RGI_OUTPUT_DIR = os.path.join(PROJECT_ROOT, config["amr_screening_dir"])

# -------------------------------------------------------------------
# CARD DB setup
# -------------------------------------------------------------------

def validate_card_db(path):
    """Ensure CARD DB has required files."""
    def has_fasta(p):
        return (
            os.path.exists(os.path.join(p, "card_reference.fasta")) or
            any(name.endswith(".fasta") for name in os.listdir(p))
        )
    if not (os.path.exists(os.path.join(path, "card.json")) and has_fasta(path)):
        sys.exit(
            f"\nERROR: CARD DB not properly prepared at {path} "
            "(missing card.json or .fasta). See pipeline README.\n"
        )

env_rgi_card = os.getenv("RGI_CARD", "").strip()
if env_rgi_card:
    RGI_CARD = env_rgi_card
elif "card_latest" in config:
    cfg_card = config["card_latest"]
    RGI_CARD = os.path.join(PROJECT_ROOT, cfg_card) if not os.path.isabs(cfg_card) else cfg_card
else:
    raise ValueError("You must set RGI_CARD in .env or card_latest in config.yaml!")

if not os.path.exists(RGI_CARD):
    sys.exit(f"\nERROR: RGI_CARD path not found: {RGI_CARD}\n")

validate_card_db(RGI_CARD)

rule rgi_reload_database:
    input:
        card_json = f"{RGI_CARD}/card.json",
        card_fasta = f"{RGI_CARD}/card_reference.fasta"
    output:
        reload_marker = f"{LOG_DIR}/rgi_reload_db.done"
    log:
        f"{LOG_DIR}/rgi/rgi_reload_db.log"
    conda:
        "../envs/rgi.yaml"
    shell:
        r"""
        set -e

        LOCALDB={RGI_CARD}

        # Standard "is file, not dir" check
        if [ -f "$LOCALDB" ]; then
            echo "ERROR: localDB exists as a file, but should be a directory." >&2
            exit 1
        fi
        mkdir -p "$LOCALDB"

        # Debug output
        echo "Current files in $LOCALDB before loading:" >> {log}
        ls -lh "$LOCALDB" >> {log}

        # Only reload if DB is not present (or is empty)
        if [ -f "$LOCALDB/card.json" ] && [ -s "$LOCALDB/card.json" ] && \
        [ -f "$LOCALDB/card_reference.fasta" ] && [ -s "$LOCALDB/card_reference.fasta" ]; then
            echo "CARD local database already present in $LOCALDB, skipping reload." >> {log}
        else
            echo "Reloading CARD into $LOCALDB..." >> {log}
            export RGI_DATA_PATH="$LOCALDB"
            rgi load --card_json {input.card_json} --card_annotation {input.card_fasta} --local >> {log} 2>&1
        fi

        touch {output.reload_marker}
        """
rule symlink_rgi_card:
    input:
        reload_marker = f"{LOG_DIR}/rgi_reload_db.done"
    output:
        done = f"{LOG_DIR}/rgi_symlink.done"
    log:
        f"{LOG_DIR}/rgi/symlink_rgi_card.log"
    shell:
        r"""
        set -o pipefail

        mkdir -p $(dirname {log})
        echo "Running on: $(hostname)" > {log}

        # Create symlink in working directory to shared CARD DB
        mkdir -p localDB
        ln -sf "{RGI_CARD}/card_reference.fasta" localDB/
        ln -sf "{RGI_CARD}/card.json" localDB/
        ln -sf "{RGI_CARD}/loaded_databases.json" localDB/
        ln -sf "{RGI_CARD}/README" localDB/
        ln -sfn "{RGI_CARD}/bwt" localDB/
        touch {output.done}
        """
rule rgi_bwt:
    input:
        R1 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz",
        R2 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz",
        symlink_marker = f"{LOG_DIR}/rgi_symlink.done"
    output:
        json = temp(f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/{{sample}}_paired.allele_mapping_data.json"),
        bai = temp(f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/{{sample}}_paired.sorted.length_100.bam.bai"),
        bam = temp(f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/{{sample}}_paired.sorted.length_100.bam"),
        allele = f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/{{sample}}_paired.allele_mapping_data.txt"
    params:
        outprefix = lambda wc: f"{CARD_RGI_OUTPUT_DIR}/{wc.sample}/{wc.sample}_paired"
    log:
        f"{LOG_DIR}/rgi/bwt_{{sample}}.log"
    threads: config.get("rgi_bwt", {}).get("threads", 4)
    conda:
        "../envs/rgi.yaml"
    shell:
        r"""
        set -o pipefail

        mkdir -p "$(dirname {log})"
        mkdir -p "$(dirname {params.outprefix})"

        # Ensure RGI_DATA_PATH is set to the localDB directory
        export RGI_DATA_PATH=localDB
        echo "RGI_DATA_PATH: $RGI_DATA_PATH" >> {log}

        rgi bwt \
            -1 {input.R1} \
            -2 {input.R2} \
            -a kma \
            -n {threads} \
            -o {params.outprefix} \
            --local \
            --clean \
            2>&1 | egrep -av '\[W::sam_parse1\]|mapped query cannot have zero coordinate' | awk 'NF' >> {log}

        rgistatus=${{PIPESTATUS[0]}}
        if [ $rgistatus -ne 0 ]; then
            echo "ERROR: rgi bwt failed with exit code $rgistatus" >> {log}
            exit $rgistatus
        fi
        """