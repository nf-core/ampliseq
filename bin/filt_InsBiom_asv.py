#!/usr/bin/env python3
##### Program description #######
#
# Title: Subsetting Insect Biome ASVs 
#
# Author(s): Lokeshwaran Manoharan
#
#
#
# Description: The ASVs that are supposed to be of particular length (418 ± nx3 and between 403 and 418), possibly does not contain any stop codon in the right reading frame. 
#  
# List of subroutines: 
#
#
#
# Overall procedure: using functions and dictionaries in python
#
# Usage: filt_InsBiom_asv.py <Seq_file> <ASV_table> <output_basename>
#
##################################

import re
import sys
import copy
import argparse

usage = '''This program takes the output from the AmpliSeq dada2 pipeline and filteres the COI ASVs and only keep the ASVs of particular length (418 ± nx3 and between 403 and 418) and does not contain any stop codon in the right reading frame'''

parser = argparse.ArgumentParser(description=usage)

parser.add_argument(
	'-f', '--fasta',
	dest='fasta',
	type=argparse.FileType('r'),
	help='Fasta file of the ASVs from the AmpliSeq pipeline',
	required=True
	)

parser.add_argument(
	'-t', '--count-table',
	dest='count',
	type=argparse.FileType('r'),
	help='Count table file of the ASVs from the AmpliSeq pipeline',
	required=True
	)

parser.add_argument(
	'-p', '--out-prefix',
	dest='prefix',
	type=str,
	help='Prefix of the output files with filtered ASVs and the corresponding count table. The output files will be in *<prefix>_filtered.* format',
	required=True
	)


args = parser.parse_args()


p1 = re.compile('\t')
p2 = re.compile('>')

def check_asv(seq):
	sub_seq = seq[1:]
	stop_list = ['TAA', 'TAG']
	sub_list = []
	for x in range(0, len(sub_seq), 3):
		sub_list.append(sub_seq[x:x+3])
	if not 402 <= len(sub_seq) <= 417:
		return False
	elif not len(sub_seq)%3 == 0:
		return False
	elif any(x in sub_list for x in stop_list):
		return False
	else:
		return True

Out_Seq = open(args.prefix +'_filtered.fna', 'w')
Out_list = open(args.prefix+'_filtered.list', 'w')
Out_table = open(args.prefix+'_filtered.table.tsv', 'w')

count_dict = {}

count = 0
for line in args.count:
	line = line.rstrip('\n')
	if count == 0:
		print(line, file = Out_table)
		count += 1
	else:
		tmp_list = re.split(p1, line)
		count_dict[tmp_list[0]] = line

for line in args.fasta:
	line = line.rstrip('\n')
	if re.match(p2, line) is not None:
		bin_head = re.sub(p2, '', line)
	elif check_asv(line):
		print('>', bin_head,'\n', line, file = Out_Seq, sep='')
		print(bin_head, file = Out_list)
		print(count_dict[bin_head], file = Out_table)

args.count.close()
args.fasta.close()
Out_Seq.close()
Out_list.close()
Out_table.close()

