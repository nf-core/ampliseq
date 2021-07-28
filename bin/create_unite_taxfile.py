#!/usr/bin/env python3
#
# Convert fasta file from UNITE into
# (i) fasta file with only sequence names and sequences
# (ii) Taxonomy file with 7 levels D_0__ to D_6__
#
# By Jeanette Tångrot 2020-09-02

#--- Import libraries, do initializations  ---#
import sys
from Bio import SeqIO

usage = """create_unite_taxfile.py <unite.fa> <unite_fasta.fa> <unite_tax.txt>
    <unite.fa> : Input. Fasta file from UNITE containing taxonomies in description line
    <unite_fasta.fa> : Output. Name of fasta file.
    <unite_tax.txt> : Output. Name of text file with taxonomies.
"""

#--- Check and read arguments ---#
if len(sys.argv) != 4:
    exit("Usage: " + usage )

fasta_in = sys.argv[1]
fasta_out = sys.argv[2]
tax_out = sys.argv[3]

#--- Read sequence file and create new records ---#
replace_dict = {';p__': ';D_1__', ';c__': ';D_2__', ';o__': ';D_3__', ';f__': ';D_4__', ';g__': ';D_5__', ';s__': ';D_6__'}

fh_fasta = open( fasta_out, mode = 'w' )
fh_tax = open( tax_out, mode = 'w' )

for entry in SeqIO.parse( fasta_in, "fasta" ):
    (name, tax) = entry.id.split('|k__')
    tax = 'D_0__' + tax
    tax = tax.replace('unidentified','')
    for n1, n2 in replace_dict.items():
        tax = tax.replace( n1, n2 )
    tax = tax.replace('|SH','_SH')
    fh_fasta.write('>' + name + '\n' + str(entry.seq).upper() + '\n')
    fh_tax.write( name + '\t' + tax + '\n' )

fh_fasta.close()
fh_tax.close()

