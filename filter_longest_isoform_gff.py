#!/usr/bin/env python3
"""Select the longest protein isoform per gene using a GFF3 annotation."""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


def parse_attributes(text):
    attributes = {}
    for item in text.rstrip().split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attributes[unquote(key)] = value
    return attributes


def attribute_values(attributes, key):
    value = attributes.get(key, "")
    return [unquote(item) for item in value.split(",") if item]


def protein_id_from_cds(attributes):
    if attributes.get("protein_id"):
        return unquote(attributes["protein_id"])

    for dbxref in attribute_values(attributes, "Dbxref"):
        if dbxref.startswith("NCBI_GP:"):
            return dbxref.removeprefix("NCBI_GP:")

    feature_id = unquote(attributes.get("ID", ""))
    if feature_id.startswith("cds-"):
        return feature_id.removeprefix("cds-")
    return None


def ncbi_geneid_from_attributes(attributes):
    """Extract the numeric NCBI GeneID from a Dbxref=GeneID:xxxx entry."""
    for dbxref in attribute_values(attributes, "Dbxref"):
        if dbxref.startswith("GeneID:"):
            return dbxref.removeprefix("GeneID:")
    return None


def load_gff_mappings(gff_path):
    transcript_to_gene = {}
    gene_names = {}
    gene_ncbi_ids = {}
    protein_to_transcripts = defaultdict(set)
    missing_protein_id = 0
    missing_parent = 0

    with gff_path.open() as gff:
        for line_number, line in enumerate(gff, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{gff_path}:{line_number}: expected 9 GFF3 columns")

            feature = fields[2]
            attributes = parse_attributes(fields[8])

            if feature == "gene":
                gene_id = unquote(attributes.get("ID", ""))
                if gene_id:
                    gene_names[gene_id] = unquote(
                        attributes.get("locus_tag", attributes.get("Name", gene_id))
                    )
                    ncbi_geneid = ncbi_geneid_from_attributes(attributes)
                    if ncbi_geneid:
                        gene_ncbi_ids[gene_id] = ncbi_geneid

            if feature in ("mRNA", "transcript") or attributes.get("gbkey") == "mRNA":
                transcript_id = unquote(attributes.get("ID", ""))
                parents = attribute_values(attributes, "Parent")
                if transcript_id and parents:
                    if len(set(parents)) > 1:
                        raise ValueError(
                            f"Transcript {transcript_id} has multiple gene parents"
                        )
                    gene_id = parents[0]
                    previous = transcript_to_gene.get(transcript_id)
                    if previous is not None and previous != gene_id:
                        raise ValueError(
                            f"Transcript {transcript_id} has multiple gene parents"
                        )
                    transcript_to_gene[transcript_id] = gene_id

            if feature != "CDS":
                continue

            protein_id = protein_id_from_cds(attributes)
            parents = attribute_values(attributes, "Parent")
            if not protein_id:
                missing_protein_id += 1
            if not parents:
                missing_parent += 1
            if protein_id and parents:
                protein_to_transcripts[protein_id].update(parents)

    protein_to_genes = defaultdict(set)
    protein_to_gene_transcripts = defaultdict(lambda: defaultdict(set))
    unresolved_transcripts = set()
    for protein_id, transcript_ids in protein_to_transcripts.items():
        for transcript_id in transcript_ids:
            gene_id = transcript_to_gene.get(transcript_id)
            if gene_id is None and transcript_id in gene_names:
                gene_id = transcript_id
            if gene_id is None:
                unresolved_transcripts.add(transcript_id)
            else:
                protein_to_genes[protein_id].add(gene_id)
                protein_to_gene_transcripts[protein_id][gene_id].add(transcript_id)

    multi_gene_proteins = {
        protein_id: sorted(gene_ids)
        for protein_id, gene_ids in protein_to_genes.items()
        if len(gene_ids) > 1
    }
    return (
        dict(protein_to_genes),
        {
            protein_id: {
                gene_id: sorted(transcript_ids)
                for gene_id, transcript_ids in genes.items()
            }
            for protein_id, genes in protein_to_gene_transcripts.items()
        },
        gene_names,
        gene_ncbi_ids,
        missing_protein_id,
        missing_parent,
        unresolved_transcripts,
        multi_gene_proteins,
    )


def iter_fasta_records(fasta_path):
    with fasta_path.open() as fasta:
        header = None
        sequence_lines = []

        def record():
            fields = header[1:].split(None, 1)
            if not fields:
                raise ValueError(f"{fasta_path}: empty FASTA identifier")
            identifier = fields[0]
            length = sum(len(line.strip()) for line in sequence_lines)
            return identifier, header, length, sequence_lines

        for line_number, line in enumerate(fasta, start=1):
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield record()
                header = line
                sequence_lines = []
            elif header is None:
                if line.strip():
                    raise ValueError(
                        f"{fasta_path}:{line_number}: sequence data before first header"
                    )
            else:
                sequence_lines.append(line)

        if header is not None:
            yield record()


def write_mapping(mapping_path, protein_to_gene, protein_to_transcripts,
                  gene_names, gene_ncbi_ids, protein_lengths, selected_assignments):
    with mapping_path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "gene_id",
                "ncbi_geneid",
                "gene_name",
                "transcript_ids",
                "protein_id",
                "protein_length",
                "selected_longest",
            )
        )
        for protein_id in sorted(protein_to_gene):
            for gene_id in sorted(protein_to_gene[protein_id]):
                writer.writerow(
                    (
                        gene_id,
                        gene_ncbi_ids.get(gene_id, ""),
                        gene_names.get(gene_id, gene_id),
                        ",".join(protein_to_transcripts[protein_id][gene_id]),
                        protein_id,
                        protein_lengths.get(protein_id, ""),
                        int((gene_id, protein_id) in selected_assignments),
                    )
                )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Select the longest protein sequence per gene by following "
            "GFF3 CDS -> transcript -> gene relationships. Output FASTA "
            "headers are the NCBI GeneID (falls back to the GFF gene ID "
            "if no GeneID Dbxref is present)."
        ),
        epilog=(
            "Example: python3 filter_longest_isoform_gff.py genomic.gff "
            "protein.faa longest_isoforms.faa --mapping gene_protein_mapping.tsv"
        ),
    )
    parser.add_argument("gff", type=Path, help="GFF3 annotation file")
    parser.add_argument("proteins", type=Path, help="protein FASTA file")
    parser.add_argument("output", type=Path, help="output FASTA file")
    parser.add_argument(
        "--mapping",
        type=Path,
        help="optional TSV mapping of genes, transcripts, proteins, and selections",
    )
    args = parser.parse_args()

    for input_path in (args.gff, args.proteins):
        if not input_path.is_file():
            parser.error(f"input file not found: {input_path}")

    output_paths = [args.output]
    if args.mapping:
        output_paths.append(args.mapping)
    protected_paths = {args.gff.resolve(), args.proteins.resolve()}
    if any(path.resolve() in protected_paths for path in output_paths):
        parser.error("output paths must not overwrite an input file")
    if len({path.resolve() for path in output_paths}) != len(output_paths):
        parser.error("FASTA output and mapping output must be different files")

    try:
        (
            protein_to_gene,
            protein_to_transcripts,
            gene_names,
            gene_ncbi_ids,
            missing_protein_id,
            missing_parent,
            unresolved_transcripts,
            multi_gene_proteins,
        ) = load_gff_mappings(args.gff)

        if not protein_to_gene:
            raise ValueError("no CDS-to-gene mappings could be read from the GFF3")

        protein_lengths = {}
        protein_sequences = {}
        best_by_gene = {}
        fasta_protein_count = 0
        for protein_id, header, length, sequence_lines in iter_fasta_records(args.proteins):
            fasta_protein_count += 1
            if protein_id in protein_lengths:
                raise ValueError(f"duplicate protein FASTA identifier: {protein_id}")
            protein_lengths[protein_id] = length
            protein_sequences[protein_id] = sequence_lines
            gene_ids = protein_to_gene.get(protein_id, set())
            for gene_id in gene_ids:
                previous = best_by_gene.get(gene_id)
                if previous is None or length > previous[0]:
                    best_by_gene[gene_id] = (length, protein_id)

        if not fasta_protein_count:
            raise ValueError(f"no FASTA records found in {args.proteins}")
        if not best_by_gene:
            raise ValueError("no protein FASTA records matched the GFF3 protein IDs")

        selected_assignments = {
            (gene_id, protein_id) for gene_id, (_, protein_id) in best_by_gene.items()
        }
        selected_proteins = {protein_id for _, protein_id in selected_assignments}

        # Translate each selected protein's header to its NCBI GeneID
        # (falls back to the GFF gene ID if no GeneID Dbxref is present).
        missing_geneid_genes = set()
        protein_to_header_geneid = {}
        for gene_id, protein_id in selected_assignments:
            ncbi_geneid = gene_ncbi_ids.get(gene_id)
            if not ncbi_geneid:
                missing_geneid_genes.add(gene_id)
                ncbi_geneid = gene_id
            protein_to_header_geneid[protein_id] = ncbi_geneid

        written_geneids = set()
        with args.output.open("w") as output:
            for protein_id, header, _, sequence_lines in iter_fasta_records(args.proteins):
                if protein_id not in selected_proteins:
                    continue
                ncbi_geneid = protein_to_header_geneid[protein_id]
                if ncbi_geneid in written_geneids:
                    continue
                written_geneids.add(ncbi_geneid)
                output.write(f">{ncbi_geneid}\n")
                for sequence_line in sequence_lines:
                    output.write(sequence_line + "\n")

        if args.mapping:
            write_mapping(
                args.mapping,
                protein_to_gene,
                protein_to_transcripts,
                gene_names,
                gene_ncbi_ids,
                protein_lengths,
                selected_assignments,
            )
    except (OSError, ValueError) as error:
        parser.error(str(error))

    matched_proteins = sum(
        protein_id in protein_lengths for protein_id in protein_to_gene
    )
    print(f"Protein FASTA records:        {fasta_protein_count:,}")
    print(f"GFF protein-to-gene mappings: {sum(map(len, protein_to_gene.values())):,}")
    print(f"Proteins matched:             {matched_proteins:,}")
    print(f"Genes retained:               {len(best_by_gene):,}")
    print(f"Proteins written:             {len(selected_proteins):,}")
    print(f"Output:                       {args.output}")
    if args.mapping:
        print(f"Mapping:                      {args.mapping}")
    if missing_protein_id:
        print(
            f"WARNING: {missing_protein_id:,} GFF CDS features had no protein ID",
            file=sys.stderr,
        )
    if missing_parent:
        print(
            f"WARNING: {missing_parent:,} GFF CDS features had no Parent",
            file=sys.stderr,
        )
    if unresolved_transcripts:
        print(
            f"WARNING: {len(unresolved_transcripts):,} CDS Parent IDs had no "
            "transcript-to-gene mapping",
            file=sys.stderr,
        )
    if missing_geneid_genes:
        print(
            f"WARNING: {len(missing_geneid_genes):,} genes had no NCBI GeneID "
            "Dbxref; used the GFF gene ID as the FASTA header instead",
            file=sys.stderr,
        )
    for protein_id, gene_ids in sorted(multi_gene_proteins.items()):
        print(
            f"WARNING: protein {protein_id} maps to multiple genes "
            f"({', '.join(gene_ids)}); considering it for each gene",
            file=sys.stderr,
        )
    if matched_proteins < len(protein_to_gene):
        print(
            f"WARNING: {len(protein_to_gene) - matched_proteins:,} GFF proteins "
            "were absent from the protein FASTA",
            file=sys.stderr,
        )
    if fasta_protein_count > matched_proteins:
        print(
            f"WARNING: {fasta_protein_count - matched_proteins:,} protein FASTA "
            "records had no GFF gene mapping",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
