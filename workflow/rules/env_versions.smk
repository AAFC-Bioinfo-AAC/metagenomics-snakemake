rule software_report:
    output:
        summary = f"{SOFTWARE_VERSIONS}/software_versions_summary.txt"
    params:
        conda_prefix = CONDA_PREFIX,
        software_versions = SOFTWARE_VERSIONS
    shell:
        r"""
        mkdir -p {params.software_versions}

        echo "==========================================================" > {output.summary}
        echo "       Conda Environment Exact Package Versions" >> {output.summary}
        echo "==========================================================" >> {output.summary}
        echo "" >> {output.summary}

        ENV_BASE="{params.conda_prefix}"

        if [ ! -d "$ENV_BASE" ]; then
            echo "ERROR: Conda env directory not found at $ENV_BASE" >> {output.summary}
            echo "Check config.yaml -> conda_prefix setting." >> {output.summary}
            exit 0
        fi

        echo "Detected Conda environments under: $ENV_BASE" >> {output.summary}
        echo "" >> {output.summary}

        for env in "$ENV_BASE"/*; do
            if [ -d "$env" ]; then
                envname=$(basename "$env")
                echo "### Environment: $envname" >> {output.summary}
                echo "Path: $env" >> {output.summary}
                echo "----------------------------------------------------------" >> {output.summary}
                conda list --prefix "$env" >> {output.summary} 2>/dev/null || echo "Could not list packages" >> {output.summary}
                echo "" >> {output.summary}
                echo "----------------------------------------------------------" >> {output.summary}
                echo "" >> {output.summary}
            fi
        done
        """
rule filter_key_bioinformatics_versions:
    input:
        f"{SOFTWARE_VERSIONS}/software_versions_summary.txt"
    output:
        f"{SOFTWARE_VERSIONS}/key_bioinformatics_software.txt"
    shell:
        r"""
        mkdir -p {SOFTWARE_VERSIONS}

        echo "==========================================================" > {output}
        echo "      Key Bioinformatics Software (Exact Versions Used)" >> {output}
        echo "==========================================================" >> {output}
        echo "" >> {output}

        KEY_TOOLS="bedtools|bowtie2|fastp|featureCounts|kraken2|megahit|pigz|prodigal|rgi|rnaQUAST|samtools|sortmerna|spades|subread"

        awk '
          /^### Environment:/ {{env=$0; print "\n" env >> "{output}"; next}}
          $1 ~ /^('"$KEY_TOOLS"')$/ {{print $0 >> "{output}"}}
        ' {input}
        """
rule filter_key_bioinformatics_html:
    input:
        f"{SOFTWARE_VERSIONS}/key_bioinformatics_software.txt"
    output:
        report(
            f"{SOFTWARE_VERSIONS}/key_bioinformatics_software.html"
        )
    shell:
        r"""
        mkdir -p {SOFTWARE_VERSIONS}
        echo "<html><body><pre>" > {output}
        cat {input} >> {output}
        echo "</pre></body></html>" >> {output}
        """

