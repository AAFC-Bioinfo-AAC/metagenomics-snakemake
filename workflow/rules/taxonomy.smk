'''
    Filename: taxonomy.smk
    Author: Katherine James-Gzyl and Devin Holman
    Date created: 2026/09/11
    Snakemake version: 9.20.0
    Python version: 3.9

'''

import os


def taxonomy_db_file(filename):
    return os.path.join(os.fspath(TAXONOMY_DB), filename)


def kraken_confidence():
    try:
        value = float(config.get("kraken2", {}).get("conf_threshold", 0.5))
    except (TypeError, ValueError) as exc:
        raise ValueError("kraken2.conf_threshold must be numeric.") from exc

    if not 0.0 <= value <= 1.0:
        raise ValueError("kraken2.conf_threshold must be between 0 and 1.")
    return value


def bracken_read_length():
    try:
        value = int(config.get("bracken", {}).get("readlen", 150))
    except (TypeError, ValueError) as exc:
        raise ValueError("bracken.readlen must be a positive integer.") from exc

    if value <= 0:
        raise ValueError("bracken.readlen must be a positive integer.")
    return value


def bracken_threshold(name, default):
    try:
        value = int(config.get("bracken", {}).get(name, default))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"bracken.{name} must be a non-negative integer.") from exc

    if value < 0:
        raise ValueError(f"bracken.{name} must be a non-negative integer.")
    return value

def bracken_domain_level():
    value = str(config.get("bracken", {}).get("domain_level", "D")).upper()
    if value not in {"D", "F1", "R1"}:
        raise ValueError(
            'bracken.domain_level must be one of "D", "F1" or "R1".'
        )
    return value


def bracken_distribution_file(wildcards):
    return taxonomy_db_file(
        f"database{bracken_read_length()}mers.kmer_distrib"
    )


rule kraken2:
    wildcard_constraints:
        sample = "[^/]+"
    input:
        hash = taxonomy_db_file("hash.k2d"),
        opts = taxonomy_db_file("opts.k2d"),
        taxo = taxonomy_db_file("taxo.k2d"),
        R1 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz",
        R2 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz"
    output:
        report = f"{KRAKEN_OUTPUT_DIR}/{{sample}}.report.txt",
        kraken = f"{KRAKEN_OUTPUT_DIR}/{{sample}}.kraken"
    log:
        f"{LOG_DIR}/kraken2/{{sample}}.log"
    conda:
        "../envs/kraken2.yaml"
    threads:
        config.get("kraken2", {}).get("threads", 2)
    params:
        db = os.fspath(TAXONOMY_DB),
        conf_threshold = lambda wildcards: kraken_confidence()
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" "$(dirname {output.report:q})"
        : > {log:q}

        kraken2 --use-names \
            --gzip-compressed \
            --threads {threads} \
            --db {params.db:q} \
            --confidence {params.conf_threshold} \
            --report-zero-counts \
            --paired {input.R1:q} {input.R2:q} \
            --report {output.report:q} \
            --output {output.kraken:q} \
            >> {log:q} 2>&1
        """


rule bracken:
    wildcard_constraints:
        sample = "[^/]+"
    input:
        distribution = bracken_distribution_file,
        report = f"{KRAKEN_OUTPUT_DIR}/{{sample}}.report.txt"
    output:
        species = f"{BRACKEN_OUTPUT_DIR}/species/{{sample}}_bracken.species.report.txt",
        genus = f"{BRACKEN_OUTPUT_DIR}/genus/{{sample}}_bracken.genus.report.txt",
        phylum = f"{BRACKEN_OUTPUT_DIR}/phylum/{{sample}}_bracken.phylum.report.txt",
        domain = f"{BRACKEN_OUTPUT_DIR}/domain/{{sample}}_bracken.domain.report.txt"
    log:
        f"{LOG_DIR}/bracken/{{sample}}.log"
    conda:
        "../envs/kraken2.yaml"
    params:
        db = os.fspath(TAXONOMY_DB),
        readlen = lambda wildcards: bracken_read_length(),
        threshold_species = lambda wildcards: bracken_threshold(
            "threshold_species", 10
        ),
        threshold_genus = lambda wildcards: bracken_threshold(
            "threshold_genus", 10
        ),
        threshold_phylum = lambda wildcards: bracken_threshold(
            "threshold_phylum", 10
        ),
        threshold_domain = lambda wildcards: bracken_threshold(
            "threshold_domain", 0
        ),
        domain_level = lambda wildcards: bracken_domain_level()
    shell:
        r"""
        set -euo pipefail

        mkdir -p \
            "$(dirname {log:q})" \
            "$(dirname {output.species:q})" \
            "$(dirname {output.genus:q})" \
            "$(dirname {output.phylum:q})" \
            "$(dirname {output.domain:q})"
        : > {log:q}

        bracken \
            -r {params.readlen} \
            -t {params.threshold_species} \
            -d {params.db:q} \
            -i {input.report:q} \
            -l S \
            -o {output.species:q} \
            >> {log:q} 2>&1

        bracken \
            -r {params.readlen} \
            -t {params.threshold_genus} \
            -d {params.db:q} \
            -i {input.report:q} \
            -l G \
            -o {output.genus:q} \
            >> {log:q} 2>&1

        bracken \
            -r {params.readlen} \
            -t {params.threshold_phylum} \
            -d {params.db:q} \
            -i {input.report:q} \
            -l P \
            -o {output.phylum:q} \
            >> {log:q} 2>&1

        bracken \
            -r {params.readlen} \
            -t {params.threshold_domain} \
            -d {params.db:q} \
            -i {input.report:q} \
            -l {params.domain_level} \
            -o {output.domain:q} \
            >> {log:q} 2>&1
        """


rule combine_bracken_outputs:
    input:
        species = expand(
            f"{BRACKEN_OUTPUT_DIR}/species/{{sample}}_bracken.species.report.txt",
            sample=SAMPLES
        ),
        genus = expand(
            f"{BRACKEN_OUTPUT_DIR}/genus/{{sample}}_bracken.genus.report.txt",
            sample=SAMPLES
        ),
        phylum = expand(
            f"{BRACKEN_OUTPUT_DIR}/phylum/{{sample}}_bracken.phylum.report.txt",
            sample=SAMPLES
        ),
        domain = expand(
            f"{BRACKEN_OUTPUT_DIR}/domain/{{sample}}_bracken.domain.report.txt",
            sample=SAMPLES
        )
    output:
        species = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_species.txt",
        genus = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_genus.txt",
        phylum = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_phylum.txt",
        domain = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_domain.txt"
    log:
        f"{LOG_DIR}/bracken/combine_bracken_outputs.log"
    conda:
        "../envs/kraken2.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p \
            "$(dirname {log:q})" \
            "$(dirname {output.species:q})"
        : > {log:q}

        echo "Combining species Bracken outputs..." >> {log:q}
        combine_bracken_outputs.py \
            --files {input.species:q} \
            --output {output.species:q} \
            >> {log:q} 2>&1

        echo "Combining genus Bracken outputs..." >> {log:q}
        combine_bracken_outputs.py \
            --files {input.genus:q} \
            --output {output.genus:q} \
            >> {log:q} 2>&1

        echo "Combining phylum Bracken outputs..." >> {log:q}
        combine_bracken_outputs.py \
            --files {input.phylum:q} \
            --output {output.phylum:q} \
            >> {log:q} 2>&1

        echo "Combining domain Bracken outputs..." >> {log:q}
        combine_bracken_outputs.py \
            --files {input.domain:q} \
            --output {output.domain:q} \
            >> {log:q} 2>&1
        """


rule clean_host_bracken:
    input:
        species = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_species.txt",
        genus = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_genus.txt",
        phylum = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_phylum.txt",
        domain = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_domain.txt"
    output:
        species = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_species_cleaned.txt",
        genus = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_genus_cleaned.txt",
        phylum = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_phylum_cleaned.txt",
        domain = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_domain_cleaned.txt"
    log:
        f"{LOG_DIR}/bracken/clean_host_bracken.log"
    conda:
        "../envs/kraken2.yaml"
    script:
        "../scripts/clean_bracken_batch.py"


rule bracken_recompute_fractions:
    input:
        species = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_species_cleaned.txt",
        genus = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_genus_cleaned.txt",
        phylum = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_phylum_cleaned.txt",
        domain = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_domain_cleaned.txt"
    output:
        species_adjusted = f"{BRACKEN_OUTPUT_DIR}/bracken_cleaned_adjusted_species.txt",
        genus_adjusted = f"{BRACKEN_OUTPUT_DIR}/bracken_cleaned_adjusted_genus.txt",
        phylum_adjusted = f"{BRACKEN_OUTPUT_DIR}/bracken_cleaned_adjusted_phylum.txt"
    log:
        f"{LOG_DIR}/bracken/recompute_fractions.log"
    conda:
        "../envs/kraken2.yaml"
    script:
        "../scripts/recompute_bracken_fractions.py"


rule bracken_extract:
    input:
        species_adjusted = f"{BRACKEN_OUTPUT_DIR}/bracken_cleaned_adjusted_species.txt",
        genus_adjusted = f"{BRACKEN_OUTPUT_DIR}/bracken_cleaned_adjusted_genus.txt",
        phylum_adjusted = f"{BRACKEN_OUTPUT_DIR}/bracken_cleaned_adjusted_phylum.txt"
    output:
        species_raw = f"{BRACKEN_OUTPUT_DIR}/bracken_species_raw_abundance.csv",
        species_rel = f"{BRACKEN_OUTPUT_DIR}/bracken_species_rel_abundance_default.csv",
        species_rel_recalc = f"{BRACKEN_OUTPUT_DIR}/bracken_species_rel_abundance_adjusted.csv",
        genus_raw = f"{BRACKEN_OUTPUT_DIR}/bracken_genus_raw_abundance.csv",
        genus_rel = f"{BRACKEN_OUTPUT_DIR}/bracken_genus_rel_abundance_default.csv",
        genus_rel_recalc = f"{BRACKEN_OUTPUT_DIR}/bracken_genus_rel_abundance_adjusted.csv",
        phylum_raw = f"{BRACKEN_OUTPUT_DIR}/bracken_phylum_raw_abundance.csv",
        phylum_rel = f"{BRACKEN_OUTPUT_DIR}/bracken_phylum_rel_abundance_default.csv",
        phylum_rel_recalc = f"{BRACKEN_OUTPUT_DIR}/bracken_phylum_rel_abundance_adjusted.csv"
    conda:
        "../envs/kraken2.yaml"
    script:
        "../scripts/extract_bracken_columns.py"
