#!/usr/bin/env python3
##### Program description #######
#
# Title: Filtering ASVs with stop codons
#
# Author(s): Lokeshwaran Manoharan
#
#
#
# Description: The ASVs are filtered if they contain stop codons in the reading frame specified.
#
# List of subroutines:
#
#
#
# Overall procedure: using functions and dictionaries in python
#
# Usage: filt_codons.py [-h] -f <Seq_file> -t <ASV_table> [-s BEGIN] [-e END] -p <output_basename>
#
##################################

import re
import sys
import copy
import argparse

usage = """This program filters ASVs that contains any stop codon or if the length from the position is not a multiple of 3 in the right reading frame."""

parser = argparse.ArgumentParser(description=usage)

parser.add_argument(
    "-f",
    "--fasta",
    dest="fasta",
    type=argparse.FileType("r"),
    help="Fasta file of the ASVs from the AmpliSeq pipeline",
    required=True,
)

parser.add_argument(
    "-t",
    "--count-table",
    dest="count",
    type=argparse.FileType("r"),
    help="Count table file of the ASVs from the AmpliSeq pipeline",
    required=False,
)

parser.add_argument(
    "-s",
    "--start-position",
    dest="begin",
    type=int,
    help="Specific position of the ASV where to start checking for codons. By default it starts from position 1",
    required=False,
    default=0,
)

parser.add_argument(
    "-e",
    "--end-position",
    dest="end",
    type=int,
    help="Specific position of the ASV where to stop checking for codons. By default it checks until the end of the sequence",
    required=False,
)

parser.add_argument(
    "-x",
    "--stop-codons",
    dest="stopcodon",
    type=str,
    help='Specific stop codons to look for. Specify multiple codons with in a comma separated list like: "TAA,TAG". By default TAA and TAG are being looked for.',
    required=False,
    default="TAA,TAG",
)

parser.add_argument(
    "-p",
    "--out-prefix",
    dest="prefix",
    type=str,
    help="Prefix of the output files with filtered ASVs and the corresponding count table. The output files will be in *<prefix>_filtered.* format",
    required=True,
)


args = parser.parse_args()

if args.begin is not None:
    begin = args.begin
else:
    begin = 1

if args.stopcodon is not None:
    stopcodon = args.stopcodon
else:
    stopcodon = "TAA,TAG"

stop_list = [item.strip() for item in stopcodon.split(",")]


def check_asv(seq, start, stop):
    sub_seq = seq[start - 1 : stop]
    sub_list = []
    for x in range(0, len(sub_seq), 3):
        sub_list.append(sub_seq[x : x + 3])
    if not len(sub_seq) % 3 == 0:
        return False
    elif any(x in sub_list for x in stop_list):
        return False
    else:
        return True


Out_Seq = open(args.prefix + "_filtered.fna", "w")
Out_list = open(args.prefix + "_filtered.list", "w")
if args.count is not None:
    Out_table = open(args.prefix + "_filtered.table.tsv", "w")
else:
    Out_table = open("empty_" + args.prefix + "_filtered.table.tsv", "w")

count_dict = {}
p1 = re.compile("\t")
p2 = re.compile(">")

if args.count is not None:
    count = 0
    for line in args.count:
        line = line.rstrip("\n")
        if count == 0:
            print(line, file=Out_table)
            count += 1
        else:
            tmp_list = re.split(p1, line)
            count_dict[tmp_list[0]] = line

for line in args.fasta:
    line = line.rstrip("\n")
    if re.match(p2, line) is not None:
        bin_head = re.sub(p2, "", line)
    else:
        if args.end is not None:
            end = args.end
        else:
            end = len(line)
        if check_asv(line, begin, end):
            print(">", bin_head, "\n", line, file=Out_Seq, sep="")
            print(bin_head, file=Out_list)
            if args.count is not None:
                print(count_dict[bin_head], file=Out_table)

if args.count is not None:
    args.count.close()
args.fasta.close()
Out_Seq.close()
Out_list.close()
Out_table.close()
