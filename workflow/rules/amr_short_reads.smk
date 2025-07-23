'''
    Filename: amr_short_reads.smk
    Author: Katherine James-Gzyl
    Date created: 2025/07/16
    Snakemake version: 9.6.0
'''
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
        R1 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
        R2 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz",
        symlink_marker = f"{LOG_DIR}/rgi_symlink.done"
    output:
        json = f"{CARD_RGI_OUTPUT_DIR}/{{sample}}_paired.allele_mapping_data.json",
        allele = f"{CARD_RGI_OUTPUT_DIR}/{{sample}}_paired.allele_mapping_data.txt"
    params:
        outprefix = lambda wc: f"{CARD_RGI_OUTPUT_DIR}/{wc.sample}_paired"
    log:
        f"{LOG_DIR}/rgi/bwt_{{sample}}.log"
    threads: 20
    conda:
        "../envs/rgi.yaml"
    shell:
        r"""
        set -o pipefail

        mkdir -p $(dirname {log})
        echo "Running on: $(hostname)" > {log}

        # Ensure RGI_DATA_PATH is set to the localDB directory
        export RGI_DATA_PATH=localDB
        echo "RGI_DATA_PATH: $RGI_DATA_PATH" >> {log}

        mkdir -p {CARD_RGI_OUTPUT_DIR}
        start_time=$(date +%s)
        echo "Started at: $(date)" >> {log}
        echo 'RGI version:' >> {log}
        rgi main --version >> {log}
        echo 'Running rgi bwt...' >> {log}

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

        end_time=$(date +%s)
        runtime=$((end_time - start_time))

        # Wall time logging
        end_time=$(date +%s)
        end_hr=$(date)
        runtime=$((end_time - start_time))
        hours=$((runtime / 3600))
        mins=$(((runtime % 3600) / 60))
        secs=$((runtime % 60))
        echo "Finished at: $end_hr" >> {log}
        echo "Wall time: $hours"h" $mins"m" $secs"s" (total $runtime seconds)" >> {log}
        """