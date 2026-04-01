# Path to the conda env YAML files, relative to the Snakefile location
ENVS_DIR = os.path.join(os.path.dirname(workflow.snakefile), "envs")
ENVS_DIR = "/fs/hpci-nas4-sci/data/aafc/projects/J-003165_abcc_rcba/workspace/metaG-testing-students/issue8_alex/metagenomics-snakemake/workflow/envs"

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

rule collect_key_software_from_yaml:
    """
    Parses all workflow/envs/*.yaml files and extracts pinned software versions.
    Only lines with a pinned version (package=X.Y.Z) under 'dependencies:' are included.
    Python-only or unpinned entries (e.g. 'python >=3.10', 'pip') are skipped.
    Output is sorted alphabetically by software name.
    """
    input:
        yamls = expand(f"{ENVS_DIR}/{{env}}.yaml", env=[
            os.path.splitext(os.path.basename(f))[0]
            for f in sorted(os.listdir(ENVS_DIR))
            if f.endswith(".yaml") or f.endswith(".yml")
        ])
    output:
        txt  = f"{SOFTWARE_VERSIONS}/key_bioinformatics_software.txt",
        html = report(f"{SOFTWARE_VERSIONS}/key_bioinformatics_software.html",
                      caption="../report/key_software.rst",
                      category="Software")
    run:
        import yaml, re

        seen = {}  # name -> (version, source_yaml)

        for yaml_path in input.yamls:
            with open(yaml_path) as fh:
                data = yaml.safe_load(fh)

            deps = data.get("dependencies", [])
            for dep in deps:
                # Skip pip sub-dict and bare strings without pinned version
                if not isinstance(dep, str):
                    continue
                # Match entries like:  packagename=1.2.3  or  packagename=1.2.3.4
                m = re.match(r'^([A-Za-z0-9_\-\.]+)=(\d[\w\.\-]*)$', dep.strip())
                if not m:
                    continue
                name, version = m.group(1).lower(), m.group(2)
                # Keep first occurrence (first yaml wins; sorted alphabetically by yaml name)
                if name not in seen:
                    seen[name] = (version, os.path.basename(yaml_path))

        rows = sorted(seen.items())  # sorted by software name

        # ── Plain text output ──────────────────────────────────────────
        with open(output.txt, "w") as out:
            out.write("=" * 58 + "\n")
            out.write("   Key Bioinformatics Software (Pinned in env YAMLs)\n")
            out.write("=" * 58 + "\n\n")
            out.write(f"{'Software':<30} {'Version':<20} {'Source YAML'}\n")
            out.write("-" * 58 + "\n")
            for name, (version, src) in rows:
                out.write(f"{name:<30} {version:<20} {src}\n")
            out.write("\n")

        # ── HTML output ────────────────────────────────────────────────
        with open(output.html, "w") as out:
            out.write("""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>Key Bioinformatics Software</title>
  <style>
    body  { font-family: Arial, sans-serif; margin: 2em; }
    h1    { color: #2c5f2e; }
    table { border-collapse: collapse; width: 60%; }
    th    { background: #2c5f2e; color: white; padding: 8px 14px; text-align: left; }
    td    { padding: 6px 14px; border-bottom: 1px solid #ddd; }
    tr:nth-child(even) { background: #f2f2f2; }
  </style>
</head>
<body>
  <h1>Key Bioinformatics Software</h1>
  <p>Pinned versions sourced from <code>workflow/envs/*.yaml</code></p>
  <table>
    <tr><th>Software</th><th>Version</th><th>Source YAML</th></tr>
""")
            for name, (version, src) in rows:
                out.write(f"    <tr><td>{name}</td><td>{version}</td><td>{src}</td></tr>\n")
            out.write("  </table>\n</body>\n</html>\n")