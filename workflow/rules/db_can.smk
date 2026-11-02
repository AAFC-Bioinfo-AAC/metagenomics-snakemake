'''
    Filename: db_can.smk
    Author: Katherine James-Gzyl and Devin Holman
    Date created: 2026/09/11
    Snakemake version: 9.20.0

'''

import os


# All samples use one DNA co-assembly and one set of gene/CGC identifiers.
DBCAN_SHARED = os.path.join(SAMPLE_DBCAN, "coassembly")
DBCAN_REFERENCE = os.path.join(DBCAN_COASSEMBLY, "final.contigs.fa")
DBCAN_ANNOTATION = os.path.join(DBCAN_SHARED, "annotation")
DBCAN_BWA_PREFIX = os.path.join(DBCAN_SHARED, "bwa_index", "coassembly")
DBCAN_BWA_INDEX = multiext(DBCAN_BWA_PREFIX, ".amb", ".ann", ".bwt", ".pac", ".sa")


def dbcan_coassembly_samples():
    excluded = config.get("dbcan_coassembly", {}).get("exclude_samples", [])
    if not isinstance(excluded, list) or any(not isinstance(x, str) for x in excluded):
        raise ValueError("dbcan_coassembly.exclude_samples must be a list of sample IDs.")
    unknown = set(excluded) - set(SAMPLE_NAMES)
    if unknown:
        raise ValueError("Unknown dbCAN co-assembly sample IDs: " + ", ".join(sorted(unknown)))
    retained = [sample for sample in SAMPLE_NAMES if sample not in set(excluded)]
    if not retained:
        raise ValueError("At least one sample must contribute reads to the dbCAN co-assembly.")
    return retained


def dbcan_coassembly_memory(wildcards, resources):
    fraction = float(config.get("dbcan_coassembly", {}).get("memory_fraction", 0.85))
    if not 0 < fraction <= 1:
        raise ValueError("dbcan_coassembly.memory_fraction must be greater than 0 and at most 1.")
    # MEGAHIT -m > 1 is an absolute limit in bytes. A fraction alone would
    # incorrectly refer to all RAM on the compute node rather than this job.
    return int(resources.mem_mb * 1_000_000 * fraction)


def bwa_mapping_threads():
    value = config.get("bwa_mem_mapping", {}).get("threads", 12)
    if isinstance(value, bool) or not isinstance(value, int) or value < 2:
        raise ValueError("bwa_mem_mapping.threads must be an integer of at least 2.")
    return value


DBCAN_COASSEMBLY_SAMPLES = dbcan_coassembly_samples()


rule dbcan_all:
    input:
        DBCAN_REFERENCE,
        f"{DBCAN_COASSEMBLY}/coassembly_samples.tsv",
        expand(
            f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/fam_abund.out",
            sample=SAMPLE_NAMES
        )


# -------------------------------------------------------------------
# Shared metagenomic co-assembly (independent of individual MAG assemblies)
# -------------------------------------------------------------------

rule megahit_coassembly:
    input:
        r1 = expand(f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz", sample=DBCAN_COASSEMBLY_SAMPLES),
        r2 = expand(f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz", sample=DBCAN_COASSEMBLY_SAMPLES)
    output:
        assembly = ensure(DBCAN_REFERENCE, non_empty=True),
        samples = f"{DBCAN_COASSEMBLY}/coassembly_samples.tsv"
    params:
        samples = DBCAN_COASSEMBLY_SAMPLES,
        min_contig_length = config.get("dbcan_coassembly", {}).get("min_contig_length", 1000),
        memory_bytes = dbcan_coassembly_memory
    threads:
        config.get("dbcan_coassembly", {}).get("threads", 24)
    resources:
        mem_mb = 256000
    log:
        f"{LOG_DIR}/dbcan/megahit_coassembly.log"
    conda:
        "../envs/megahit.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log:q})" "$(dirname {output.assembly:q})"
        : > {log:q}

        # Use the output filesystem for work files: large co-assemblies must
        # not depend on a small or periodically cleaned node-local /tmp.
        work=$(mktemp -d "$(dirname {output.assembly:q})/.megahit.XXXXXX")
        trap 'rm -rf -- "$work"' EXIT
        mkdir -p "$work/reads" "$work/tmp"
        r1=( {input.r1:q} )
        r2=( {input.r2:q} )
        samples=( {params.samples:q} )
        forward=()
        reverse=()
        printf 'sample\tfastq_1\tfastq_2\n' > "$work/samples.tsv"
        for i in "${{!r1[@]}}"; do
            [[ -s "${{r1[$i]}}" && -s "${{r2[$i]}}" ]] || {{
                echo "Missing or empty read file for ${{samples[$i]}}" >> {log:q}
                exit 1
            }}
            # Short symlink names also support paths containing spaces/commas.
            ln -s -- "${{r1[$i]}}" "$work/reads/${{i}}_R1.fastq.gz"
            ln -s -- "${{r2[$i]}}" "$work/reads/${{i}}_R2.fastq.gz"
            forward+=( "reads/${{i}}_R1.fastq.gz" )
            reverse+=( "reads/${{i}}_R2.fastq.gz" )
            printf '%s\t%s\t%s\n' "${{samples[$i]}}" "${{r1[$i]}}" "${{r2[$i]}}" >> "$work/samples.tsv"
        done
        r1_list=$(IFS=,; echo "${{forward[*]}}")
        r2_list=$(IFS=,; echo "${{reverse[*]}}")

        (
            cd "$work"
            megahit -1 "$r1_list" -2 "$r2_list" \
                -t {threads} -m {params.memory_bytes} \
                --min-contig-len {params.min_contig_length} \
                --out-prefix final --tmp-dir tmp -o run
        ) >> {log:q} 2>&1

        [[ -s "$work/run/final.contigs.fa" ]] || {{
            echo "ERROR: MEGAHIT produced no contigs; no placeholder assembly was created." >> {log:q}
            exit 1
        }}
        mv -- "$work/run/final.contigs.fa" {output.assembly:q}
        mv -- "$work/samples.tsv" {output.samples:q}
        """


rule bwa_index_coassembly:
    input:
        assembly = DBCAN_REFERENCE
    output:
        DBCAN_BWA_INDEX
    params:
        prefix = DBCAN_BWA_PREFIX
    threads: 1
    log:
        f"{LOG_DIR}/dbcan/bwa_index_coassembly.log"
    conda:
        "../envs/bwa.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log:q})" "$(dirname {params.prefix:q})"
        bwa index -p {params.prefix:q} {input.assembly:q} > {log:q} 2>&1
        """


# -------------------------------------------------------------------
# Gene prediction
# -------------------------------------------------------------------

rule pyrodigal:
    input:
        assembly = DBCAN_REFERENCE,
        normalize_gff = os.path.join(WORKFLOW_DIR, "scripts", "normalize_dbcan_gff.py")
    output:
        gff = f"{DBCAN_SHARED}/coassembly_genes.gff",
        faa = f"{DBCAN_SHARED}/coassembly_proteins.faa",
        nuc_file = f"{DBCAN_SHARED}/coassembly.cds"
    log:
        f"{LOG_DIR}/dbcan/prodigal/coassembly.log"
    threads:
        config.get("pyrodigal", {}).get("threads", 8)
    conda:
        "../envs/pyrodigal.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" "$(dirname {output.gff:q})"
        : > {log:q}

        raw_gff=$(mktemp "$(dirname {output.gff:q})/.genes.XXXXXX.gff")
        trap 'rm -f -- "$raw_gff"' EXIT

        pyrodigal \
            -p meta \
            -i {input.assembly:q} \
            -d {output.nuc_file:q} \
            -a {output.faa:q} \
            -f gff \
            -o "$raw_gff" \
            -j {threads} \
            >> {log:q} 2>&1

        python {input.normalize_gff:q} \
            --proteins {output.faa:q} --gff "$raw_gff" --output {output.gff:q} \
            >> {log:q} 2>&1
        [[ -s {output.faa:q} && -s {output.nuc_file:q} ]]
        """


# -------------------------------------------------------------------
# Shared run_dbCAN analyses
#
# The default workflow runs easy_substrate once. The two standalone rules
# are available as explicit alternatives and do not run by default.
# -------------------------------------------------------------------

rule cazyme_annotation:
    input:
        faa = f"{DBCAN_SHARED}/coassembly_proteins.faa",
        dbcan_db = DB_CAN_DB_PATH
    output:
        overview = f"{DBCAN_SHARED}/annotation_cazyme/overview.tsv"
    params:
        outdir = f"{DBCAN_SHARED}/annotation_cazyme"
    log:
        f"{LOG_DIR}/dbcan/cazyme_annotation/coassembly.log"
    threads:
        config.get("cazyme_annotation", {}).get("threads", 8)
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" {params.outdir:q}
        : > {log:q}

        run_dbcan CAZyme_annotation \
            --input_raw_data {input.faa:q} \
            --mode protein \
            --output_dir {params.outdir:q} \
            --db_dir {input.dbcan_db:q} \
            --threads {threads} \
            >> {log:q} 2>&1
        """


rule cgc_calling:
    input:
        faa = f"{DBCAN_SHARED}/coassembly_proteins.faa",
        gff = f"{DBCAN_SHARED}/coassembly_genes.gff",
        dbcan_db = DB_CAN_DB_PATH
    output:
        overview = f"{DBCAN_SHARED}/annotation_cgc/overview.tsv",
        cgc = f"{DBCAN_SHARED}/annotation_cgc/cgc.gff",
        cgc_table = f"{DBCAN_SHARED}/annotation_cgc/cgc_standard_out.tsv"
    params:
        outdir = f"{DBCAN_SHARED}/annotation_cgc"
    log:
        f"{LOG_DIR}/dbcan/cgc_calling/coassembly.log"
    threads:
        config.get("cgc_calling", {}).get("threads", 8)
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" {params.outdir:q}
        : > {log:q}

        run_dbcan easy_CGC \
            --input_raw_data {input.faa:q} \
            --mode protein \
            --input_gff {input.gff:q} \
            --gff_type prodigal \
            --output_dir {params.outdir:q} \
            --db_dir {input.dbcan_db:q} \
            --threads {threads} \
            >> {log:q} 2>&1
        """


rule substrate_prediction:
    input:
        faa = f"{DBCAN_SHARED}/coassembly_proteins.faa",
        gff = f"{DBCAN_SHARED}/coassembly_genes.gff",
        dbcan_db = DB_CAN_DB_PATH
    output:
        overview = f"{DBCAN_SHARED}/annotation/overview.tsv",
        cgc_table = f"{DBCAN_SHARED}/annotation/cgc_standard_out.tsv",
        substrate = f"{DBCAN_SHARED}/annotation/substrate_prediction.tsv",
        subfamily = f"{DBCAN_SHARED}/annotation/dbCANsub_hmm_results.tsv"
    params:
        outdir = f"{DBCAN_SHARED}/annotation"
    log:
        f"{LOG_DIR}/dbcan/substrate_prediction/coassembly.log"
    threads:
        config.get("substrate_prediction", {}).get("threads", 8)
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" {params.outdir:q}
        : > {log:q}

        run_dbcan easy_substrate \
            --input_raw_data {input.faa:q} \
            --mode protein \
            --input_gff {input.gff:q} \
            --gff_type prodigal \
            --output_dir {params.outdir:q} \
            --db_dir {input.dbcan_db:q} \
            --threads {threads} \
            >> {log:q} 2>&1
        """


# -------------------------------------------------------------------
# Read mapping and gene-depth calculation
# -------------------------------------------------------------------

rule bwa_mem_mapping:
    input:
        assembly = DBCAN_REFERENCE,
        index = DBCAN_BWA_INDEX,
        R1 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz",
        R2 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz"
    output:
        bam = f"{DBCAN_SHARED}/samples/{{sample}}/mapping/{{sample}}.bam",
        bai = f"{DBCAN_SHARED}/samples/{{sample}}/mapping/{{sample}}.bam.bai"
    params:
        index_prefix = DBCAN_BWA_PREFIX
    log:
        f"{LOG_DIR}/dbcan/bwa_mem_mapping/{{sample}}.log"
    threads:
        bwa_mapping_threads()
    conda:
        "../envs/bwa.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" "$(dirname {output.bam:q})"
        : > {log:q}

        if (( {threads} < 2 )); then
            echo "ERROR: bwa_mem_mapping requires at least 2 threads; received {threads}." >> {log:q}
            exit 1
        fi

        # Reserve approximately one quarter of the CPU budget for SAMtools.
        sort_threads=$(( {threads} / 4 ))
        if (( sort_threads < 1 )); then
            sort_threads=1
        fi
        bwa_threads=$(( {threads} - sort_threads ))

        # Reserve one CPU within the SAMtools budget for its main thread.
        sort_extra=$(( sort_threads - 1 ))

        echo "[$(date)] Total threads: {threads}; BWA: $bwa_threads; SAMtools budget: $sort_threads; SAMtools workers: $sort_extra" \
            >> {log:q}

        work=$(mktemp -d "$(dirname {output.bam:q})/.sort.XXXXXX")
        trap 'rm -rf -- "$work"' EXIT

        bwa mem \
            -t "$bwa_threads" \
            {params.index_prefix:q} \
            {input.R1:q} \
            {input.R2:q} \
            2>> {log:q} \
        | samtools sort \
            -@ "$sort_extra" -m 1G -T "$work/sort" \
            -o {output.bam:q} \
            - \
            >> {log:q} 2>&1

        samtools quickcheck -v {output.bam:q} >> {log:q} 2>&1
        samtools index {output.bam:q} {output.bai:q} >> {log:q} 2>&1
        """


rule dbcan_depth:
    input:
        gff = f"{DBCAN_SHARED}/coassembly_genes.gff",
        bam = f"{DBCAN_SHARED}/samples/{{sample}}/mapping/{{sample}}.bam",
        bai = f"{DBCAN_SHARED}/samples/{{sample}}/mapping/{{sample}}.bam.bai"
    output:
        depth_file = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/{{sample}}.depth.txt"
    log:
        f"{LOG_DIR}/dbcan/dbcan_depth/{{sample}}.log"
    params:
        overlap_base_ratio = config.get("dbcan_depth", {}).get(
            "overlap_base_ratio", 0.2
        ),
        mapping_quality = config.get("dbcan_depth", {}).get(
            "mapping_quality", 30
        ),
        identity = config.get("dbcan_depth", {}).get("identity", 0.98)
    threads:
        config.get("dbcan_depth", {}).get("threads", 16)
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" "$(dirname {output.depth_file:q})"
        : > {log:q}

        dbcan_utils cal_coverage \
            -g {input.gff:q} \
            -i {input.bam:q} \
            -o {output.depth_file:q} \
            -t {threads} \
            --overlap_base_ratio {params.overlap_base_ratio} \
            --mapping_quality {params.mapping_quality} \
            --identity {params.identity} \
            >> {log:q} 2>&1
        """


# -------------------------------------------------------------------
# Abundance calculation
# -------------------------------------------------------------------

rule get_abundances_rpm:
    input:
        overview = f"{DBCAN_SHARED}/annotation/overview.tsv",
        cgc_table = f"{DBCAN_SHARED}/annotation/cgc_standard_out.tsv",
        substrate = f"{DBCAN_SHARED}/annotation/substrate_prediction.tsv",
        subfamily = f"{DBCAN_SHARED}/annotation/dbCANsub_hmm_results.tsv",
        depth_file = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/{{sample}}.depth.txt"
    output:
        cazy_fam_ab = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/fam_abund.out",
        cazy_subfam_ab = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/subfam_abund.out",
        EC_number = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/EC_abund.out",
        substrate_ab = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/fam_substrate_abund.out",
        cgc_ab = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/CGC_abund.out",
        substrate_ho = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/CGC_substrate_PUL_homology.out",
        cgc_substrate_voting = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/CGC_substrate_majority_voting.out"
    params:
        dbcan_dir = f"{DBCAN_SHARED}/annotation",
        ab_dir = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund",
        marker_file = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/did_not_run_get_abundances_rpm.txt"
    log:
        f"{LOG_DIR}/dbcan/get_abundances_rpm/{{sample}}.log"
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" {params.ab_dir:q}
        : > {log:q}
        rm -f -- {params.marker_file:q} {output:q}

        # Return success only when a tabular file contains at least one
        # non-empty row after its header.
        has_data_rows() {{
            local input_file="$1"
            [[ -s "$input_file" ]] && \
                awk 'NR > 1 && NF {{ found = 1; exit }} END {{ exit !found }}' \
                    "$input_file"
        }}

        require_output() {{
            local output_file="$1"
            if [[ ! -e "$output_file" ]]; then
                echo "ERROR: Expected output was not created: $output_file" \
                    >> {log:q}
                exit 1
            fi
        }}

        record_skip() {{
            printf '%s\n' "$1" >> {params.marker_file:q}
        }}

        (
            cd {params.ab_dir:q}

            if has_data_rows {input.overview:q}; then
                dbcan_utils fam_abund \
                    -bt {input.depth_file:q} \
                    -i {params.dbcan_dir:q} \
                    -a RPM \
                    >> {log:q} 2>&1

                require_output {output.cazy_fam_ab:q}
                require_output {output.cazy_subfam_ab:q}
                require_output {output.EC_number:q}

                dbcan_utils fam_substrate_abund \
                    -bt {input.depth_file:q} \
                    -i {params.dbcan_dir:q} \
                    -a RPM \
                    >> {log:q} 2>&1

                require_output {output.substrate_ab:q}
            else
                : > {output.cazy_fam_ab:q}
                : > {output.cazy_subfam_ab:q}
                : > {output.EC_number:q}
                : > {output.substrate_ab:q}
                record_skip \
                    "Family-level abundance calculations were skipped because overview.tsv contained no annotation rows."
            fi

            if has_data_rows {input.cgc_table:q}; then
                dbcan_utils CGC_abund \
                    -bt {input.depth_file:q} \
                    -i {params.dbcan_dir:q} \
                    -a RPM \
                    >> {log:q} 2>&1

                require_output {output.cgc_ab:q}
            else
                : > {output.cgc_ab:q}
                record_skip \
                    "CGC abundance calculation was skipped because cgc_standard_out.tsv contained no CGC rows."
            fi

            if has_data_rows {input.substrate:q}; then
                dbcan_utils CGC_substrate_abund \
                    -bt {input.depth_file:q} \
                    -i {params.dbcan_dir:q} \
                    -a RPM \
                    >> {log:q} 2>&1

                require_output {output.substrate_ho:q}
                require_output {output.cgc_substrate_voting:q}
            else
                : > {output.substrate_ho:q}
                : > {output.cgc_substrate_voting:q}
                record_skip \
                    "CGC substrate abundance calculations were skipped because substrate_prediction.tsv contained no prediction rows."
            fi
        )

        """
