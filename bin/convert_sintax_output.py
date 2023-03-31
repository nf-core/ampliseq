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
# Usage: convert_sintax_output.py -i <sintax_results.tsv> -o <output_filename> -t <taxlevels>
##################################

import re
import argparse

usage = """This program takes the output of SINTAX and creates taxonomy file in a format suitable for the ampliseq pipeline"""

parser = argparse.ArgumentParser(description=usage)

parser.add_argument(
    "-i", "--infile", dest="infile", metavar="INFILE", type=argparse.FileType("r"), help="SINTAX output", required=True
)

parser.add_argument(
    "-o",
    "--outfile",
    dest="outfile",
    metavar="OUTFILE",
    type=argparse.FileType("w"),
    help="Final Taxonomy in tsv format",
)

parser.add_argument(
    "-t",
    "--taxlevels",
    dest="taxlevels",
    metavar="TAXLEVELS",
    help="Taxonomy levels to use, comma separated list",
    default="",
)

args = parser.parse_args()


def complete_list(some_list, target_len):
    return some_list[:target_len] + ["NA"] * (target_len - len(some_list))


header = ["ASV"] + args.taxlevels.split(",")
num_taxa = len(header) - 1
print("\t".join(header), file=args.outfile)

for line in args.infile:
    line = line.rstrip("\n")
    tmp_list = line.split("\t")
    if tmp_list[3] != "":
        annot = re.sub(".\:", "", tmp_list[3])
        tmp_list1 = annot.split(",")
        if len(tmp_list1) != num_taxa:
            tmp_list1 = complete_list(tmp_list1, num_taxa)
            print(tmp_list[0], "\t".join(tmp_list1), sep="\t", file=args.outfile)
        else:
            print(tmp_list[0], "\t".join(tmp_list1), sep="\t", file=args.outfile)
    else:
        print(tmp_list[0], "\t".join(["NA"] * num_taxa), sep="\t", file=args.outfile)

args.outfile.close()
args.infile.close()
