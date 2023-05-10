#!/usr/bin/env python3
##### Program description #######
#
# Title: Convert sintax output to DADA2-like tsv format
#
# Author(s): Lokeshwaran Manoharan, Jeanette Tångrot
#
# Description:
#
# List of subroutines:
#
# Usage: convert_sintax_output.py -i <sintax_results.tsv> -f <fastafile> -o <output_filename> -t <taxlevels> -d <db_version>
##################################

import re
import argparse

usage = """This program takes the output of SINTAX and creates taxonomy file in a format suitable for the ampliseq pipeline"""

parser = argparse.ArgumentParser(description=usage)

parser.add_argument(
    "-i",
    "--infile",
    dest="infile",
    metavar="INFILE",
    type=argparse.FileType("r"),
    help="SINTAX output",
    required=True,
)

parser.add_argument(
    "-f",
    "--fasta",
    dest="fastafile",
    metavar="FASTAFILE",
    type=argparse.FileType("r"),
    help="FASTA file with sequences",
    required=True,
)

parser.add_argument(
    "-o",
    "--outfile",
    dest="outfile",
    metavar="OUTFILE",
    type=argparse.FileType("w"),
    help="Final Taxonomy in tsv format",
    default="sintax_taxonomy.tsv",
)

parser.add_argument(
    "-t",
    "--taxlevels",
    dest="taxlevels",
    metavar="TAXLEVELS",
    help="Taxonomy levels to use, comma separated list",
    default="",
)

parser.add_argument(
    "-d",
    "--dbversion",
    dest="db",
    metavar="DBVERSION",
    help="Taxonomy database name",
    default="",
)

args = parser.parse_args()


def complete_list(some_list, target_len):
    return some_list[:target_len] + [""] * (target_len - len(some_list))


# Find out if taxonomies contain SH groups
dbtype = "default"
if "UNITE" in args.db:
    dbtype = "UNITE"

# Read fasta file and store as dictionary
# If multiple sequences have the same name only the first will be stored
seqs = dict()
seq = ""
name = ""
for line in args.fastafile:
    if line.startswith(">"):
        if seq != "" and name != "":
            if name not in seqs:
                seqs[name] = seq
        seq = ""
        name = line.lstrip(">").rstrip()
    else:
        seq = seq + line.rstrip("\n")
if seq != "" and name != "" and name not in seqs:
    seqs[name] = seq

# Print header to outfile
if args.taxlevels != "":
    if dbtype == "UNITE":
        header = ["ASV_ID"] + args.taxlevels.split(",") + ["SH", "confidence", "sequence"]
    else:
        header = ["ASV_ID"] + args.taxlevels.split(",") + ["confidence", "sequence"]
else:
    header = ["ASV_ID", "confidence", "sequence"]
num_taxa = len(header) - 3
num_cols = num_taxa
if dbtype == "UNITE":
    num_taxa = num_taxa - 1
print("\t".join(header), file=args.outfile)

# Read sintax file, parse results, and write taxonomies together with sequence to outfile
for line in args.infile:
    line = line.rstrip("\n")
    tmp_list = line.split("\t")
    sequence = seqs[tmp_list[0]] if tmp_list[0] in seqs else ""
    if tmp_list[3] != "":
        annot = re.sub("[dkpcofgs]\:", "", tmp_list[3])
        annot_list = annot.split(",")
        raw_annot = re.sub("[dkpcofgs]\:", "", tmp_list[1])
        raw_annot_list = raw_annot.split(",")
        confidence = raw_annot_list[len(annot_list) - 1].split("(")[1]
        confidence = confidence[:-1]
        if len(annot_list) != num_taxa:
            annot_list = complete_list(annot_list, num_cols)
            print(tmp_list[0], "\t".join(annot_list), confidence, sequence, sep="\t", file=args.outfile)
        else:
            if dbtype == "UNITE":
                (species, sh) = annot_list[num_taxa - 1].rsplit("_", 1)
                annot_list[num_taxa - 1] = species
                annot_list.append(sh)
            print(tmp_list[0], "\t".join(annot_list), confidence, sequence, sep="\t", file=args.outfile)
    else:
        print(tmp_list[0], "\t".join([""] * num_cols), "", sequence, sep="\t", file=args.outfile)


args.outfile.close()
args.infile.close()
