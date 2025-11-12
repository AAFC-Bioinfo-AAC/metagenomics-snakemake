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

The **Metagenomics Snakemake pipeline** was designed for paired-end reads from host-associated metagenomic samples that were sequenced on an Illumina platform. The reproducible workflow is modular, allowing modules to be omitted provided cleaned are provided as input. Utilizing the **Snakemake** workflow, the pipeline automates read quality control, filtering, and host/Phix decontamination, as well as downstream modules that obtain biological information from the cleaned reads. The modules provide taxonomic abundance tables, antibiotic resistance profiles, higher order functional pathways, and metagenome-assembled genome (MAG)s assembled from individual samples.

The pipeline consists of these modules:

**Pre-processing**: *fastp* is used for quality control, trimming, and filtering paired-end reads. Host and PhiX decontamination is performed using *Bowtie2* and *Samtools* by aligning the trimmed reads to the host and PhiX reference genome and retaining only those reads that do not align. For low-biomass samples, two iterations of the host removal may be required. Inclution of the host genome into the taxonomy database enables assesment of host contamination.

**Taxonomy profiling**: *Kraken2* and *Bracken* are used to generate abundance tables classifying the short reads at species, genus, and phylum level across all samples. The *Kraken2*-formatted database includes host taxonomy, and reads classified as host are removed prior to calculating the relative abundance. Instructions for modifying host parameters are provided in the user guide.

**Antimicrobial resistance profiling**: *RGI* (Resistance Gene Identifier) is used to map the short reads to *CARD* (Comprehensive Antibiotic Resistance Database). This produces a report containing the potential antibiotic resistance genes for each sample.

**Functional pathways**: High-level functional pathways are generated using the *KEGG* (Kyoto Encyclopedia of Genes and Genomes) database. The cleaned paired-end reads are concatenated and aligned to the *KEGG* protein database using *DIAMOND*. Resulting alignments are used to produce a KO (KEGG Orthology) table, which with some formatting, is used as input for *MinPath*. *MinPath* infers the minimum set of pathways required to describe the genes present. Using the *MinPath*-confirmed pathways, the pathways are summarized into high-level categories with the KEGG BRITE hierarchy. For each sample, the abundance of each higher-level pathway is expressed as reads-per-kilobase (RPK) and copies per million (CPM).

**MAGs from individual samples**: Each sample is assembled using *MEGAHIT*. The resulting assemblies are indexed, and the reads are mapped back with *Bowtie2* and *Samtools*. *MetaBAT* bins the assembled contigs into putative genomes, which are evaluated for completeness and contamination using *CheckM2*.

---

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
