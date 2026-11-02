"""Match Pyrodigal GFF gene IDs to its protein FASTA IDs for dbCAN counting.

Prodigal-compatible GFFs may use IDs such as 1_1 or FASTA IDs such as
k141_42_1. dbcan_utils joins counts to annotations by ID, so they must agree.
Coordinates and all other GFF attributes are preserved.
"""

import argparse
import re


def normalize(proteins, gff, output):
    identifiers = {}
    protein_ids = set()
    with open(proteins, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            protein_id = line[1:].split()[0]
            match = re.search(r"(?:^|[;\s])ID=([^;\s]+)", line)
            if not match:
                raise ValueError(f"No Pyrodigal ID attribute in protein header: {line.strip()}")
            original_id = match.group(1)
            if original_id in identifiers or protein_id in protein_ids:
                raise ValueError(f"Duplicate Pyrodigal/protein ID: {original_id} / {protein_id}")
            identifiers[original_id] = protein_id
            protein_ids.add(protein_id)
    if not identifiers:
        raise ValueError("Pyrodigal produced no protein sequences.")

    seen = set()
    with open(gff, encoding="utf-8") as source, open(output, "w", encoding="utf-8") as target:
        for line in source:
            if line.startswith("#") or not line.strip():
                target.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                raise ValueError("Expected nine-column Pyrodigal CDS records.")
            attrs = fields[8].split(";")
            ids = [attr[3:] for attr in attrs if attr.startswith("ID=")]
            if len(ids) != 1:
                raise ValueError(f"Missing, duplicate or unmatched GFF ID: {ids}")
            original_id = ids[0]
            protein_id = original_id if original_id in protein_ids else identifiers.get(original_id)
            if protein_id is None or protein_id in seen:
                raise ValueError(f"Missing, duplicate or unmatched GFF ID: {ids}")
            seen.add(protein_id)
            fields[8] = ";".join(
                "ID=" + protein_id if attr.startswith("ID=") else attr
                for attr in attrs
            )
            target.write("\t".join(fields) + "\n")
    if seen != protein_ids:
        raise ValueError("Some predicted proteins have no corresponding GFF gene.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    normalize(args.proteins, args.gff, args.output)
