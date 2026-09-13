'''
    Filename: preprocessing.smk
    Author: Katherine James-Gzyl and Devin Holman
    Date created: 2026/09/11
    Snakemake version: 9.20.0
'''

rule fastp_pe:
    input:
        fastq1 = lambda wc: SAMPLES[wc.sample]["fastq_1"],
        fastq2 = lambda wc: SAMPLES[wc.sample]["fastq_2"]
    output:
        r1   = temp(f"{TRIMMED_DIR}/{{sample}}_r1.fastq.gz"),
        r2   = temp(f"{TRIMMED_DIR}/{{sample}}_r2.fastq.gz"),
        u1   = temp(f"{TRIMMED_DIR}/{{sample}}_u1.fastq.gz"),
        u2   = temp(f"{TRIMMED_DIR}/{{sample}}_u2.fastq.gz"),
        html = temp(f"{TRIMMED_DIR}/{{sample}}.fastp.html"),
        json = temp(f"{TRIMMED_DIR}/{{sample}}.fastp.json")
    log:
        f"{LOG_DIR}/fastp/{{sample}}.fastp.log"
    params:
        cut_tail = "--cut_tail" if config.get("fastp", {}).get("cut_tail", True) else "",
        cut_front = "--cut_front" if config.get("fastp", {}).get("cut_front", True) else "",
        detect_adapter = "--detect_adapter_for_pe" if config.get("fastp", {}).get("detect_adapter_for_pe", True) else "",
        cut_mean_quality = config.get("fastp", {}).get("cut_mean_quality", 20),
        cut_window_size = config.get("fastp", {}).get("cut_window_size", 4),
        qualified_quality_phred = config.get("fastp", {}).get("qualified_quality_phred", 15),
        length_required = config.get("fastp", {}).get("length_required", 100)
    threads: config.get("fastp", {}).get("threads", 2)
    conda:
        "../envs/fastp.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.r1:q})" "$(dirname {log:q})"

        fastp \
            --in1 {input.fastq1:q} \
            --in2 {input.fastq2:q} \
            --out1 {output.r1:q} \
            --out2 {output.r2:q} \
            --unpaired1 {output.u1:q} \
            --unpaired2 {output.u2:q} \
            {params.cut_tail} \
            {params.cut_front} \
            {params.detect_adapter} \
            --cut_mean_quality {params.cut_mean_quality} \
            --cut_window_size {params.cut_window_size} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --length_required {params.length_required} \
            --json {output.json:q} \
            --html {output.html:q} \
            --thread {threads} \
            > {log:q} 2>&1
        """


rule bowtie2_align:
    input:
        r1 = f"{TRIMMED_DIR}/{{sample}}_r1.fastq.gz",
        r2 = f"{TRIMMED_DIR}/{{sample}}_r2.fastq.gz",
        idx = BOWTIE_INDEX_FILES
    output:
        bam = temp(f"{TRIMMED_DIR}/bam/{{sample}}.bam")
    log:
        f"{LOG_DIR}/bowtie2/{{sample}}.log"
    params:
        rg_id = lambda wc: wc.sample,
        rg_sm = lambda wc: f"SM:{wc.sample}"
    threads: config.get("bowtie2_align", {}).get("threads", 12)
    conda:
        "../envs/bowtie2.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.bam:q})" "$(dirname {log:q})"

        if (( {threads} < 3 )); then
            echo "ERROR: bowtie2_align requires at least 3 threads; received {threads}." >> {log:q}
            exit 1
        fi

        # bowtie2, samtools view, and samtools sort run concurrently in this pipe.
        # Reserve one CPU for samtools view and one for samtools sort's main thread.
        bt2_threads=$(( {threads} / 2 ))
        [ $bt2_threads -lt 1 ] && bt2_threads=1
        sort_extra=$(( {threads} - bt2_threads - 2 ))

        if (( sort_extra > 0 )); then
            sort_threads="-@ $sort_extra"
        else
            sort_threads=""
        fi

        bowtie2 \
            -x {BOWTIE_INDEX:q} \
            -1 {input.r1:q} \
            -2 {input.r2:q} \
            --threads "$bt2_threads" \
            --rg-id {params.rg_id:q} \
            --rg {params.rg_sm:q} \
            2>> {log:q} \
        | samtools view -u 2>> {log:q} \
        | samtools sort $sort_threads -o {output.bam:q} - 2>> {log:q}
        """


rule extract_unmapped_fastq:
    input:
        bam = f"{TRIMMED_DIR}/bam/{{sample}}.bam"
    output:
        r1 = protected(f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz"),
        r2 = protected(f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz")
    log:
        f"{LOG_DIR}/bedtools/{{sample}}.log"
    threads: config.get("extract_unmapped_fastq", {}).get("threads", 5)
    conda:
        "../envs/bedtools.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.r1:q})" "$(dirname {log:q})"

        if (( {threads} < 5 )); then
            echo "ERROR: extract_unmapped_fastq requires at least 5 threads; received {threads}." >> {log:q}
            exit 1
        fi

        # samtools view, samtools sort, bedtools and two compressors
        # run concurrently. Reserve one CPU for each.
        sort_extra=$(( {threads} - 5 ))
        if (( sort_extra > 0 )); then
            sort_threads="-@ $sort_extra"
        else
            sort_threads=""
        fi

        tmpbase="${{TMPDIR:-/tmp}}"
        job_id="${{SLURM_JOB_ID:-manual}}"
        mkdir -p "$tmpbase"
        tmpjob=$(mktemp -d "${{tmpbase%/}}/bam2fq_${{job_id}}_XXXXXX")

        pid1=""
        pid2=""

        cleanup() {{
            for pid in "$pid1" "$pid2"; do
                if [[ -n "$pid" ]]; then
                    kill "$pid" 2>/dev/null || true
                    wait "$pid" 2>/dev/null || true
                fi
            done
            rm -rf -- "$tmpjob"
        }}
        trap cleanup EXIT

        # Explicit child PIDs let us detect a failed compressor and wait
        # for both complete gzip streams before publishing the outputs.
        mkfifo "$tmpjob/r1.fifo" "$tmpjob/r2.fifo"

        pigz -p 1 --fast \
            < "$tmpjob/r1.fifo" \
            > "$tmpjob/r1.fastq.gz" \
            2>> {log:q} &
        pid1=$!

        pigz -p 1 --fast \
            < "$tmpjob/r2.fifo" \
            > "$tmpjob/r2.fastq.gz" \
            2>> {log:q} &
        pid2=$!

        # -f 12 retains read pairs for which both the read and its mate are unmapped.
        # -F 256 excludes secondary alignments. The name sort is required by
        # bedtools bamtofastq when producing paired FASTQ files with -fq2.
        samtools view -u -f 12 -F 256 {input.bam:q} 2>> {log:q} \
        | samtools sort -n $sort_threads \
            -T "$tmpjob/{wildcards.sample}_sort_tmp" \
            -O BAM - 2>> {log:q} \
        | bedtools bamtofastq -i - 2>> {log:q} \
            -fq "$tmpjob/r1.fifo" \
            -fq2 "$tmpjob/r2.fifo"

        wait "$pid1"
        pid1=""

        wait "$pid2"
        pid2=""

        mv -- "$tmpjob/r1.fastq.gz" {output.r1:q}
        mv -- "$tmpjob/r2.fastq.gz" {output.r2:q}
        """
