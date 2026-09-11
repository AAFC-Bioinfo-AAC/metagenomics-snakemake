<!-- omit in toc -->
# METAGENOMICS SNAKEMAKE PIPELINE

[![FR](https://img.shields.io/badge/lang-FR-yellow.svg)](README_FR.md)
[![EN](https://img.shields.io/badge/lang-EN-blue.svg)](README.md)

---

<!-- omit in toc -->
## Table of Contents

- [About](#about)
- [Documentation](#documentation)
- [Acknowledgements](#acknowledgements)
- [Security](#security)
- [License](#license)

---

## About

The **Metagenomics Snakemake pipeline** is a reproducible workflow for paired-end Illumina shotgun metagenomic reads from high-biomass, host-associated samples. It automates read quality control, trimming, filtering, and removal of host and PhiX sequences, followed by several downstream analyses. The workflow rules are organized into modules, allowing users to run selected analysis targets when their required upstream files are available. Most downstream modules use host-depleted paired reads, while the CAZyme module also requires per-sample assemblies. Depending on the analyses run, the pipeline produces taxonomic abundance tables, antimicrobial resistance gene profiles, functional pathway profiles, metagenome-assembled genomes (MAGs), and carbohydrate-active enzyme (CAZyme) annotations and abundance tables.

The pipeline contains the following modules:

- **Pre-processing:** *fastp* performs read quality assessment, adapter trimming, quality trimming, and length filtering. The trimmed reads are aligned against a combined host and PhiX reference using *Bowtie2*. *SAMtools* and *BEDTools* are then used to retain paired reads for which neither mate aligns to the reference.
- **Taxonomic profiling:** *Kraken2* classifies the cleaned reads against a user-specified Kraken2-formatted reference database; *Bracken* estimates abundances at the domain, phylum, genus, and species levels. Counts assigned to configured host taxa are removed before relative abundances are recalculated. The workflow provides raw counts, default Bracken relative abundances, and adjusted relative abundances normalized to the total number of bacterial and archaeal reads. The database name, taxonomy release, construction date, and compatible Bracken read-length distribution should be documented in the configuration or user guide.
  
- **Antimicrobial resistance gene profiling:** *RGI* (Resistance Gene Identifier) uses KMA to map the cleaned paired-end reads against reference sequences from the *Comprehensive Antibiotic Resistance Database* (CARD). The per-sample reports list putative antimicrobial resistance gene matches. These sequence-based predictions do not establish gene expression or phenotypic resistance.
  
- **Functional pathways:** Cleaned R1 and R2 read files are concatenated and searched against the *Kyoto Encyclopedia of Genes and Genomes* (KEGG) protein database using *DIAMOND* `blastx`. Gene hits are assigned to KEGG Orthology (KO) identifiers and KO abundance tables are generated. *MinPath* infers a parsimonious set of pathways consistent with the detected KOs. Each inferred pathway is assigned a pathway name and corresponding top-level and sublevel categories from the KEGG BRITE hierarchy. Per-sample pathway abundances are reported as raw counts, reads per kilobase (RPK), and counts per million reads (CPM).
  
- **MAGs from individual samples:** Each sample is assembled separately using *MEGAHIT*. A checkpoint retains only assemblies that satisfy configurable thresholds for total assembled length and number of contigs. For each retained assembly, reads from the same sample are mapped back to the assembled contigs using *Bowtie2* and *SAMtools*. *MetaBAT2* groups the contigs into putative genome bins and *CheckM2* estimates the completeness and contamination of the resulting bins.
  
- **CAZyme annotation:** A checkpoint retains only samples with a non-empty per-sample assembly. *Pyrodigal* predicts genes and proteins from each retained assembly. The predicted proteins are analyzed using [run_dbCAN](https://run-dbcan.readthedocs.io/en/latest/index.html) to identify CAZymes, CAZyme gene clusters (CGCs), and predicted substrates. Cleaned reads are mapped back to the corresponding assembly using *BWA-MEM* and *SAMtools*. *dbcan_utils* then generates gene-depth files and calculates CAZyme abundances as reads per million (RPM).
<br>
Reference databases are supplied separately and are not distributed with the pipeline. For reproducibility, users should record the database names, releases, construction dates, and relevant compatibility settings, including the Bracken read-length distribution.

## Documentation

For technical details, including installation and usage instructions, please see the [**`User Guide`**](./docs/user-guide.md).

---

## Acknowledgements

- **Credits**: This project was developed at the *Lacombe Research and Development Centre, Agriculture & Agri-Food Canada (AAFC)* by **Katherine James-Gzyl** and assisted by **Devin Holman** and **Arun Kommadath**.

- **Citation**: To cite this project, click the **`Cite this repository`** button on the right-hand sidebar

- **Contributing**: Contributions are welcome! Please review the guidelines in [CONTRIBUTING.md](CONTRIBUTING.md) and ensure you adhere to our [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) to foster a respectful and inclusive environment.

- **References**: For a list of key resources used here, see [REFERENCES.md](REFERENCES.md)

---

## Security  

⚠️ Do not post any security issues on the public repository! Please report them as described in [SECURITY.md](SECURITY.md)

---

## License

See the [LICENSE](LICENSE) file for details. Visit [LicenseHub](https://licensehub.org) or [tl;drLegal](https://www.tldrlegal.com/) to view a plain-language summary of this license.

**Copyright ©** His Majesty the King in Right of Canada, as represented by the Minister of Agriculture and Agri-Food, 2025.

---
