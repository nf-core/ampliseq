#!/usr/bin/env python3

import argparse
import re


def read_fasta_sequences(handle):
    seqs = {}
    seq_id = ""
    chunks = []
    for line in handle:
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            if seq_id and seq_id not in seqs:
                seqs[seq_id] = "".join(chunks)
            seq_id = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.strip())
    if seq_id and seq_id not in seqs:
        seqs[seq_id] = "".join(chunks)
    return seqs


def parse_lca_tokens(raw_taxonomy):
    if not raw_taxonomy:
        return []
    return [tok.strip() for tok in re.split(r"[;,]", raw_taxonomy) if tok.strip()]


def rank_key_to_name(token_key):
    key = token_key.strip().lower()
    rank_map = {
        "d": "Kingdom",
        "k": "Kingdom",
        "kingdom": "Kingdom",
        "p": "Phylum",
        "phylum": "Phylum",
        "c": "Class",
        "class": "Class",
        "o": "Order",
        "order": "Order",
        "f": "Family",
        "family": "Family",
        "g": "Genus",
        "genus": "Genus",
        "s": "Species",
        "species": "Species",
    }
    return rank_map.get(key, "")


def parse_taxonomy(raw_taxonomy, ordered_taxlevels):
    tax = {level: "" for level in ordered_taxlevels}
    tokens = parse_lca_tokens(raw_taxonomy)

    sequential_idx = 0
    for token in tokens:
        token = re.sub(r"\([^)]*\)$", "", token).strip()
        if not token:
            continue

        if "__" in token:
            key, value = token.split("__", 1)
            rank = rank_key_to_name(key)
        elif ":" in token:
            key, value = token.split(":", 1)
            rank = rank_key_to_name(key)
        else:
            # Fallback: fill taxlevels left to right when no prefixes are present
            if sequential_idx >= len(ordered_taxlevels):
                continue
            rank = ordered_taxlevels[sequential_idx]
            value = token
            sequential_idx += 1

        value = value.strip()
        if rank in tax and value:
            tax[rank] = value

    return tax


def infer_taxonomy_field(fields):
    if len(fields) < 2:
        return ""
    for field in fields[1:]:
        if "__" in field or ":" in field or ";" in field or "," in field:
            return field
    return fields[-1]


def main():
    parser = argparse.ArgumentParser(description="Convert VSEARCH lcaout to a DADA2-like taxonomy TSV.")
    parser.add_argument("-i", "--infile", required=True, type=argparse.FileType("r"))
    parser.add_argument("-f", "--fasta", required=True, type=argparse.FileType("r"))
    parser.add_argument("-o", "--outfile", required=True, type=argparse.FileType("w"))
    parser.add_argument(
        "-t",
        "--taxlevels",
        default="Kingdom,Phylum,Class,Order,Family,Genus,Species",
        help="Comma separated list of taxonomy levels to emit",
    )
    args = parser.parse_args()

    taxlevels = [x.strip() for x in args.taxlevels.split(",") if x.strip()]
    seqs = read_fasta_sequences(args.fasta)

    header = ["ASV_ID"] + taxlevels + ["confidence", "sequence"]
    print("\t".join(header), file=args.outfile)

    for line in args.infile:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        asv_id = fields[0].strip()
        raw_taxonomy = infer_taxonomy_field(fields).strip()
        taxonomy = parse_taxonomy(raw_taxonomy, taxlevels)
        confidence = ""
        sequence = seqs.get(asv_id, "")
        values = [asv_id] + [taxonomy[level] for level in taxlevels] + [confidence, sequence]
        print("\t".join(values), file=args.outfile)

    args.infile.close()
    args.fasta.close()
    args.outfile.close()


if __name__ == "__main__":
    main()
