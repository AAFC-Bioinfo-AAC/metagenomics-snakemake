DBCAN_ENVS = ["pyrodigal", "dbcan", "bwa"]

localrules: prewarm_dbcan_env, prewarm_dbcan_gate

rule prewarm_dbcan_env:
    output:
        f"{LOG_DIR}/envs/dbcan/{{env}}.prewarmed"
    conda:
        "../envs/{env}.yaml"
    shell:
        "mkdir -p {LOG_DIR}/envs/dbcan && touch {output}"
rule prewarm_dbcan_gate:
    input:
        expand(f"{LOG_DIR}/envs/dbcan/{{env}}.prewarmed", env=DBCAN_ENVS)
    output:
        touch(f"{LOG_DIR}/envs/conda_gate_dbcan.txt")
    shell:
        "mkdir -p {LOG_DIR}/envs && touch {output}"
checkpoint nonempty_assemblies:
    input:
        assemblies = expand(f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa", sample=SAMPLE_NAMES),
        gate = f"{LOG_DIR}/envs/conda_gate_dbcan.txt"
    output:
        f"{SAMPLE_ASSEMBLY}/nonempty_assemblies.txt"
    run:
        import os
        valid_samples = []
        for infile in input.assemblies:
            if os.path.isfile(infile) and os.path.getsize(infile) > 0:
                sample = os.path.basename(infile).replace("_assembly.contigs.fa", "")
                valid_samples.append(sample)
        with open(output[0], "w") as out:
            out.write("\n".join(valid_samples) + ("\n" if valid_samples else ""))
rule pyrodigal:
    input:
        assembly = f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa"
    output:
        gff = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_genes.gff", 
        faa = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_proteins.faa", 
        nuc_file = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}.cds"
    log:
        f"{LOG_DIR}/dbcan/prodigal/{{sample}}.log"
    threads: config.get("pyrodigal", {}).get("threads", 8)
    conda:
        "../envs/pyrodigal.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        mkdir -p "$(dirname {output.gff})"
        
        pyrodigal \
            -p meta \
            -i {input.assembly} \
            -d {output.nuc_file} \
            -a {output.faa} \
            -f gff \
            -o {output.gff} \
            -j {threads} \
            >> {log} 2>&1
        """
rule cazyme_annotation:#This rule might be redundant and possible could be removed. 
    input:
        faa = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_proteins.faa",
        dbcan_db = DB_CAN_DB_PATH
    output:
        out_dir = directory(f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_cazyme"),
        overview = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_cazyme/overview.tsv"
    params:
        outdir = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_cazyme"
    log:
        f"{LOG_DIR}/dbcan/cazyme_annotation/{{sample}}.log"
    threads: config.get("cazyme_annotation", {}).get("threads", 8)
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        mkdir -p {params.outdir}
        
        run_dbcan CAZyme_annotation \
            --input_raw_data {input.faa} \
            --mode protein \
            --output_dir {params.outdir} \
            --db_dir {input.dbcan_db} \
            --threads {threads} \
            >> {log} 2>&1
        """
rule cgc_calling:#This rule might be redundant and possible could be removed.
    input:
        faa = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_proteins.faa",
        gff = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_genes.gff",
        dbcan_db = DB_CAN_DB_PATH
    output:
        out_dir = directory(f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_pul"),
        cgc = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_pul/cgc.gff"
    params:
        outdir = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_pul"
    log:
        f"{LOG_DIR}/dbcan/cgc_calling/{{sample}}.log"
    threads: config.get("cgc_calling", {}).get("threads", 8)
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        mkdir -p {params.outdir}

        run_dbcan easy_CGC \
            --input_raw_data {input.faa} \
            --mode protein \
            --input_gff {input.gff} \
            --output_dir {params.outdir} \
            --db_dir {input.dbcan_db} \
            --threads {threads} \
            >> {log} 2>&1
        """
rule substrate_prediction:
    input:
        faa = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_proteins.faa",
        gff = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_genes.gff",
        dbcan_db = DB_CAN_DB_PATH
    output:
        out_dir = directory(f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_dbcan"),
        overview = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_dbcan/overview.tsv",
        substrate = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_dbcan/substrate_prediction.tsv"
    params:
        outdir = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_dbcan"
    log:
        f"{LOG_DIR}/dbcan/substrate_prediction/{{sample}}.log"
    threads: config.get("substrate_prediction", {}).get("threads", 8)
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        mkdir -p {params.outdir}

        run_dbcan easy_substrate \
            --input_raw_data {input.faa} \
            --mode protein \
            --input_gff {input.gff} \
            --output_dir {params.outdir} \
            --db_dir {input.dbcan_db} \
            --threads {threads} \
            >> {log} 2>&1
        """
rule bwa_mem_mapping:
    input:
        assembly = f"{SAMPLE_ASSEMBLY}/{{sample}}_assembly.contigs.fa",
        R1 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz",
        R2 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz"
    output:
        bam = f"{SAMPLE_DBCAN}/{{sample}}/mapping/{{sample}}.bam",
        bai = f"{SAMPLE_DBCAN}/{{sample}}/mapping/{{sample}}.bam.bai"
    log:
        f"{LOG_DIR}/dbcan/bwa_mem_mapping/{{sample}}.log"
    threads: config.get("bwa_mem_mapping", {}).get("threads", 12) #thread splitting
    shadow: "minimal"
    conda:
        "../envs/bwa.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        mkdir -p "$(dirname {output.bam})"
        
        # Calculate thread distribution: ~25% to samtools, rest to BWA
        sort_threads=$(( {threads}/4 ))              # 25% for samtools (floor)
        bwa_threads=$(( {threads} - sort_threads ))  # BWA gets the remainder
        
        # Ensure minimum of 1 thread each
        [ $bwa_threads -lt 1 ] && bwa_threads=1
        [ $sort_threads -lt 1 ] && sort_threads=1
        
        echo "[$(date)] Total threads: {threads}, BWA: $bwa_threads, samtools: $sort_threads" >> {log}
        
        bwa index {input.assembly} 2>> {log}
        
        bwa mem -t $bwa_threads {input.assembly} {input.R1} {input.R2} 2>> {log} | \
            samtools sort -@ $sort_threads -o {output.bam} - >> {log} 2>&1
        
        samtools index {output.bam} >> {log} 2>&1
        """
rule dbcan_depth:
    input:
        gff = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_genes.gff",
        bam = f"{SAMPLE_DBCAN}/{{sample}}/mapping/{{sample}}.bam"
    output:
        depth_file = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/{{sample}}.depth.txt"
    log:
        f"{LOG_DIR}/dbcan/dbcan_depth/{{sample}}.log"
    params:
        overlap_base_ratio = config.get("dbcan_depth", {}).get("overlap_base_ratio", 0.2), 
        mapping_quality = config.get("dbcan_depth", {}).get("mapping_quality", 30), 
        identity = config.get("dbcan_depth", {}).get("identity", 0.98) 
    conda:
        "../envs/dbcan.yaml"
    threads: config.get("dbcan_depth",{}).get("threads",16)
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        mkdir -p "$(dirname {output.depth_file})"

        dbcan_utils cal_coverage \
            -g {input.gff} \
            -i {input.bam} \
            -o {output.depth_file} \
            -t {threads} \
            --overlap_base_ratio {params.overlap_base_ratio} \
            --mapping_quality {params.mapping_quality} \
            --identity {params.identity} >> {log} 2>&1
        """
rule get_abundances_rpm:
    input:
        overview = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_dbcan/overview.tsv",
        depth_file = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/{{sample}}.depth.txt"
    output:
        cazy_fam_ab = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/fam_abund.out",
        cazy_subfam_ab = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/subfam_abund.out",
        EC_number = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/EC_abund.out",
        substrate_ab = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/fam_substrate_abund.out",
        cgc_ab = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/CGC_abund.out",
        substrate_ho = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/CGC_substrate_PUL_homology.out",
        cgc_substrate_voting = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund/CGC_substrate_majority_voting.out"
    params:
        dbcan_dir = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_dbcan",
        ab_dir = f"{SAMPLE_DBCAN}/{{sample}}/{{sample}}_abund"
    log:
        f"{LOG_DIR}/dbcan/get_abundances_rpm/{{sample}}.log"
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"

        (
          cd "{params.ab_dir}"

          # run tools (no stdout redirects)
        dbcan_utils fam_abund -bt "{input.depth_file}" -i "{params.dbcan_dir}" -a RPM >> "{log}" 2>&1
        dbcan_utils fam_substrate_abund -bt "{input.depth_file}" -i "{params.dbcan_dir}" -a RPM >> "{log}" 2>&1
        dbcan_utils CGC_abund -bt "{input.depth_file}" -i "{params.dbcan_dir}" -a RPM >> "{log}" 2>&1
        dbcan_utils CGC_substrate_abund -bt "{input.depth_file}" -i "{params.dbcan_dir}" -a RPM >> "{log}" 2>&1
        )
        """
