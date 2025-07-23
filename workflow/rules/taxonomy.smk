'''
    Filename: taxonomy.smk
    Author: Katherine James-Gzyl
    Date created: 2025/07/16
    Snakemake version: 9.6.0
'''
rule kraken2:
    wildcard_constraints:
        sample = '[^/]+'
    input:
        ref = f"{TAXONOMY_DB}",
        R1 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
        R2 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz"
    output:
        report = f"{KRAKEN_OUTPUT_DIR}/{{sample}}.report.txt",
        kraken = f"{KRAKEN_OUTPUT_DIR}/{{sample}}.kraken" 
    log:
        f"{LOG_DIR}/kraken2/{{sample}}.log"
    conda:
        "../envs/kraken2.yaml"
    threads: 2
    params:
        conf_threshold = "0.5"
    shell:
        r"""
        set -euo pipefail

        # Wall time tracking
        start_time=$(date +%s)
        start_hr=$(date)
        echo "Started at: $start_hr" > {log}

        echo 'Kraken2 version:' >> {log}
        kraken2 --version >> {log}
        # echo 'Database is GTDB 220' >> {log} #uncomment if there is a text file with the version
        # cat {input.ref}/VERSION.txt >> {log} #uncomment if there is a text file with the version 

        kraken2 --use-names \
            --threads {threads} \
            --db {input.ref} \
            --confidence {params.conf_threshold} \
            --report-zero-counts \
            --paired {input.R1} {input.R2} \
            --report {output.report} \
            --output {output.kraken} \
            &>> {log}
            
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
rule bracken:
    wildcard_constraints:
        sample = '[^/]+'
    input:
        ref = f"{TAXONOMY_DB}",
        report = f"{KRAKEN_OUTPUT_DIR}/{{sample}}.report.txt"
    output:
        species = f"{KRAKEN_OUTPUT_DIR}/species/{{sample}}_bracken.species.report.txt",
        genus = f"{KRAKEN_OUTPUT_DIR}/genus/{{sample}}_bracken.genus.report.txt",
        phylum = f"{KRAKEN_OUTPUT_DIR}/phylum/{{sample}}_bracken.phylum.report.txt"
    log:
        f"{LOG_DIR}/bracken/{{sample}}.log"
    conda:
        "../envs/kraken2.yaml"
    threads: 2
    params:
        readlen = 150
    shell:
        r"""
        set -euo pipefail

        # Wall time tracking
        start_time=$(date +%s)
        start_hr=$(date)
        echo "Started at: $start_hr" > {log}

        echo 'bracken version:' >> {log}
        bracken --version >> {log}

        mkdir -p $(dirname {output.species}) $(dirname {output.genus}) $(dirname {output.phylum})

        bracken -r {params.readlen} -t {threads} -d {input.ref} -i {input.report} -l S -o {output.species} &>> {log}

        bracken -r {params.readlen} -t {threads} -d {input.ref} -i {input.report} -l G -o {output.genus} &>> {log}

        bracken -r {params.readlen} -t {threads} -d {input.ref} -i {input.report} -l P -o {output.phylum} &>> {log}

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
rule combine_bracken_outputs:
    input:
        species = expand(f"{KRAKEN_OUTPUT_DIR}/species/{{sample}}_bracken.species.report.txt", sample=SAMPLES),
        genus = expand(f"{KRAKEN_OUTPUT_DIR}/genus/{{sample}}_bracken.genus.report.txt", sample=SAMPLES),
        phylum = expand(f"{KRAKEN_OUTPUT_DIR}/phylum/{{sample}}_bracken.phylum.report.txt", sample=SAMPLES)
    output:
        species = f"{KRAKEN_OUTPUT_DIR}/merged_abundance_species.txt",
        genus = f"{KRAKEN_OUTPUT_DIR}/merged_abundance_genus.txt",
        phylum = f"{KRAKEN_OUTPUT_DIR}/merged_abundance_phylum.txt"
    log:
        f"{LOG_DIR}/bracken/combine_bracken_outputs.log"
    conda:
        "../envs/kraken2.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Combining species bracken outputs..." > {log}
        combine_bracken_outputs.py --files {input.species} --output {output.species} &>> {log}

        echo "Combining genus bracken outputs..." >> {log}
        combine_bracken_outputs.py --files {input.genus} --output {output.genus} &>> {log}

        echo "Combining phylum bracken outputs..." >> {log}
        combine_bracken_outputs.py --files {input.phylum} --output {output.phylum} &>> {log}
        """
rule bracken_extract:
    input:
        species_table = f"{KRAKEN_OUTPUT_DIR}/merged_abundance_species.txt",
        genus_table = f"{KRAKEN_OUTPUT_DIR}/merged_abundance_genus.txt",
        phylum_table = f"{KRAKEN_OUTPUT_DIR}/merged_abundance_phylum.txt"
    output:
        species_raw = f"{KRAKEN_OUTPUT_DIR}/Bracken_species_raw_abundance.csv",
        species_rel = f"{KRAKEN_OUTPUT_DIR}/Bracken_species_relative_abundance.csv", 
        genus_raw = f"{KRAKEN_OUTPUT_DIR}/Bracken_genus_raw_abundance.csv",
        genus_rel = f"{KRAKEN_OUTPUT_DIR}/Bracken_genus_relative_abundance.csv",
        phylum_raw = f"{KRAKEN_OUTPUT_DIR}/Bracken_phylum_raw_abundance.csv",
        phylum_rel = f"{KRAKEN_OUTPUT_DIR}/Bracken_phylum_relative_abundance.csv"
    conda:
        "../envs/kraken2.yaml"
    script:
        "../scripts/extract_bracken_columns.py"