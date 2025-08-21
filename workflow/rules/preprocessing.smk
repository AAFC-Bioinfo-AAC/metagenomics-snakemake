'''
    Filename: preprocessing.smk
    Author: Katherine James-Gzyl
    Date created: 2025/07/24
    Snakemake version: 9.6.0
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
        cut_tail = "--cut_tail" if config["fastp"].get("cut_tail", True) else "",
        cut_front = "--cut_front" if config["fastp"].get("cut_front", True) else "",
        detect_adapter = "--detect_adapter_for_pe" if config["fastp"].get("detect_adapter_for_pe", True) else "",
        cut_mean_quality = config["fastp"].get("cut_mean_quality", 20),
        cut_window_size = config["fastp"].get("cut_window_size", 4),
        qualified_quality_phred = config["fastp"].get("qualified_quality_phred", 15),
        length_required = config["fastp"].get("length_required", 100)
    threads: config["fastp"].get("threads", 4)
    conda: "../envs/fastp.yaml"
    shell:
        r"""
        fastp \
            --in1 {input.fastq1} \
            --in2 {input.fastq2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --unpaired1 {output.u1} \
            --unpaired2 {output.u2} \
            {params.cut_tail} {params.cut_front} {params.detect_adapter} \
            --cut_mean_quality {params.cut_mean_quality} \
            --cut_window_size {params.cut_window_size} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --length_required {params.length_required} \
            --json {output.json} \
            --html {output.html} \
            --thread {threads} \
            > {log} 2>&1
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
        extra= lambda wc: f"-R '@RG\\tID:{wc.sample}\\tSM:{wc.sample}'"
    threads: config.get("bowtie2_align", {}).get("threads", 12)
    conda:
        "../envs/bowtie2.yaml"
    shell:
        r"""
        # Split threads safely
        bt2_threads=$(( ({threads}+1)/2 ))
        sort_threads=$(( ({threads}+2)/3 ))
        view_threads=$(( {threads} - bt2_threads - sort_threads ))
        # Ensure no thread count is 0
        [ $bt2_threads -lt 1 ] && bt2_threads=1
        [ $sort_threads -lt 1 ] && sort_threads=1
        [ $view_threads -lt 1 ] && view_threads=1

        set -euo pipefail

        bowtie2 -x {BOWTIE_INDEX} -1 {input.r1} -2 {input.r2} --threads $bt2_threads {params.extra} 2>> {log} \
        | samtools view -u -@ $view_threads 2>> {log} \
        | samtools sort -@ $sort_threads -o {output.bam} 2>> {log}
        """  
rule extract_unmapped_fastq:
    input:
        bam=f"{TRIMMED_DIR}/bam/{{sample}}.bam"
    output:
        r1=protected(f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz"),
        r2=protected(f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz")

    log:
        f"{LOG_DIR}/bedtools/{{sample}}.log"
    threads: config["extract_unmapped_fastq"].get("threads", 24)
    conda:
        "../envs/bedtools.yaml"
    shell:
        r"""
        # Split threads safely
        samtools_threads=$(( ({threads}*4 + 4)/5 ))   # ~80% for samtools
        pigz_threads=$(( {threads} - samtools_threads )) # rest for pigz
        # Ensure no thread count is 0
        [ $samtools_threads -lt 1 ] && samtools_threads=1
        [ $pigz_threads -lt 1 ] && pigz_threads=1
        
        TMPDIR="${{TMPDIR:-/tmp}}"
        JOB_ID="${{SLURM_JOB_ID:-manual}}"
        TMPJOB=$(mktemp -d "${{TMPDIR}}/bam2fq_${{JOB_ID}}_XXXXXX" 2>> {log} || echo "/tmp/bam2fq_${{JOB_ID}}_manualtmp")

        samtools view -u -f 12 -F 256 --threads $samtools_threads {input.bam} 2>> {log} \
        | samtools sort -n --threads $samtools_threads -T $TMPJOB/{wildcards.sample}_sort_tmp -O BAM - 2>> {log} \
        | bedtools bamtofastq -i - 2>> {log} \
            -fq >(pigz -p $pigz_threads --fast > {output.r1}) \
            -fq2 >(pigz -p $pigz_threads --fast > {output.r2})

        echo "Temporary directory to be removed: $TMPJOB" >> {log}
        rm -rf "$TMPJOB"
        """