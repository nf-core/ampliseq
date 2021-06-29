#!/usr/bin/env python3
#@author Jeanette Tangrot
# Takes one TSV taxonomy file from DADA2 and a sequence fasta file,
# adds sequence to taxonomy based on ASV_ID

import pandas as pd
import sys, os 

# Argument check
if len(sys.argv) != 3:
    exit("Usage: add_full_sequence_to_taxfile.py <ASV_tax.tsv> <ASV_seqs.fasta>")

# Read tsv and remove sequence column
taxfile = sys.argv[1]
tax = pd.read_csv(taxfile, sep='\t', header=0)
tax.drop(columns='sequence', inplace=True)

# Read fasta file and store as data frame
seqs = pd.DataFrame(columns=["id","sequence"])
seq = ""
name = ""
with open(sys.argv[2], 'r') as reader:
    for line in reader:
        if line.startswith('>'):
            if (seq != "" and name != ""):
                seqs = seqs.append({'id':name, 'sequence': seq}, ignore_index=True)
                seq = ""
            name = line.lstrip('>').rstrip('\s+*\n')
        else:
            seq = seq + line.rstrip('\n')
if (seq != "" and name != ""):
    seqs = seqs.append({'id':name, 'sequence': seq}, ignore_index=True)

# Create name of results file
outfile = taxfile.replace("ASV_ITS_", "ASV_")

# Join taxonomy and full sequence, write to file
tax = tax.set_index('ASV_ID').join(seqs.set_index('id'), how='outer')
tax.to_csv(outfile, sep='\t',na_rep="NA", index_label="ASV_ID")

