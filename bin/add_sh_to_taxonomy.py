#!/usr/bin/env python3
#@author Jeanette Tångrot

# Adds SH information to ASV table based on vsearch usearch_global results in blast6 format.
#
# Usage: add_sh_to_taxonomy.py <seq2sh.tsv> <SHs.tax> <tax.tsv> <blastout.tab> <outfile> [--species]
#   Input:
#          <seq2sh.tsv>   : List with sequence to SH matchings, SH level 1.5, for each sequence in
#                           the database used by vsearch. Each row should contain database sequence
#                           name, SH name, and taxon number for the SH, separated by tabs.
#          <SHs.tax>      : List of complete taxonomy for each SH.
#          <tax.tsv>      : ASV taxonomy table to add SH assignments to
#          <blastout.tab> : Results from vsearch usearch_global in blast6 format
#          <outfile>      : Name of results file, i.e. the updated ASV taxonomy table
#   Options:
#          --species : Whether or not to include species information. Default: false

import sys
import pandas as pd

with_species = False

# Argument check
if len(sys.argv) != 6 and len(sys.argv) != 7:
    exit("Usage: add_sh_to_taxonomy.py <seq2sh.tsv> <SHs.tax> <tax.tsv> <blastout.tab> <outfile> [--species]")
if len(sys.argv) == 7:
    if sys.argv[6] == '--species':
        with_species = True
    else:
        exit("Usage: add_sh_to_taxonomy.py <seq2sh.tsv> <SHs.tax> <tax.tsv> <blastout.tab> <outfile> [--species]")
outfile = sys.argv[5]

# Set number of ranks - include species if --species flag is used
if with_species:
    num_ranks = 8
    tax_entries = ['Domain','Kingdom','Phylum','Class','Order','Family','Genus','Species','SH','confidence']
else:
    num_ranks = 7
    tax_entries = ['Domain','Kingdom','Phylum','Class','Order','Family','Genus','SH','confidence']
    
# Read sequence to SH matchings
seq2sh = pd.read_csv(sys.argv[1], sep='\t', header=None, index_col=0, skiprows=None)

# Read SH taxonomies
# Columns:
# SH  taxonid  kingdom  phylum  class  order  family  genus  species
shtax = pd.read_csv(sys.argv[2], sep='\t', header=None, index_col=0, skiprows=None)
# Replace taxonid with Domain = "Eukaryota"
shtax.loc[:,1] = 'Eukaryota'

# Read taxonomy table
# ASV_ID  Domain  Kingdom Phylum  Class   Order   Family  Genus   confidence      sequence
taxtable = pd.read_csv(sys.argv[3], sep='\t', header=0)
# Add SH slot to table:
# ASV_ID  Domain  Kingdom Phylum  Class   Order   Family  Genus  SH confidence      sequence
taxtable.insert(num_ranks+1,"SH","", allow_duplicates=False)

# Go through vsearch matches and update taxonomy for those entries
fh = open( sys.argv[4], mode = 'r' )
prev_ASV = fh.readline().split()[0]
fh.seek(0)
matches = []
maxid = -1
maxlen = -1
for row in fh:
    [ASV, match, pid, alen, therest] = row.split(maxsplit=4)
    pid = float(pid)
    alen = float(alen)
    
    if ASV != prev_ASV:
        SH = ""
        tax = ""
        conf = 0.0
        for m in matches:
            matchparts = m[0].split('|')
            try:
                new_SH = seq2sh.loc[ matchparts[1] ][1]
            except KeyError:
                print( "WARNING: " + matchparts[1] + " not in seq2SH list", file=sys.stderr )
                new_SH = ""
            if ( pd.isna( new_SH ) ) :
                print( "WARNING: no SH reported for " + matchparts[1], file=sys.stderr )
                new_SH = ""                
            if SH != "" and new_SH != SH :
                SH = ""
                tax = ""
                break
            elif new_SH != "":
                SH = new_SH
                try:
                    tax = list(shtax.loc[ SH ])
                except KeyError:
                    print( "WARNING: no taxonomy found for " + SH, file=sys.stderr )
                    tax = [""]*num_ranks
                conf = m[1]/100.0
        if SH != "":
            tax_list = tax[0:num_ranks] + [SH] + [conf]
            taxtable.loc[ taxtable['ASV_ID'] == prev_ASV, tax_entries] = tax_list
        prev_ASV = ASV
        maxid = -1
        maxlen = -1
        matches = []
    if match != "*":
        if pid > maxid:
            maxid = pid
            maxlen = alen
            matches = []
            matches.append([match, pid, alen])
        elif pid == maxid and alen > maxlen:
            maxlen = pid
            matches = []
            matches.append([match, pid, alen])
        elif pid == maxid and alen == maxlen:
            matches.append([match, pid, alen])

if match != "*":        # Take care of last row/ASV in match file
    SH = ""
    tax = ""
    conf = 0.0
    for m in matches:
        matchparts = m[0].split('|')
        try:
            new_SH = seq2sh.loc[ matchparts[1] ][1]
        except KeyError:
            print( "WARNING: " + matchparts[1] + " not in seq2SH list", file=sys.stderr )
            new_SH = ""
        if ( pd.isna( new_SH ) ) :
            print( "WARNING: no SH reported for " + matchparts[1], file=sys.stderr )
            new_SH = ""                
        if SH != "" and new_SH != SH :
            SH = ""
            tax = ""
            break
        elif new_SH != "":
            SH = new_SH
            try:
                tax = list(shtax.loc[ SH ])
            except KeyError:
                print( "WARNING: no taxonomy found for " + SH, file=sys.stderr )
                tax = [""]*num_ranks
            conf = m[1]/100.0
    if SH != "":
        tax_list = tax[0:num_ranks] + [SH] + [conf]
        taxtable.loc[ taxtable['ASV_ID'] == prev_ASV, tax_entries] = tax_list

   
# Write new taxtable, with SH and new taxonomy added if found
taxtable.to_csv(outfile, sep="\t", na_rep="", float_format="%.2f", index=False)

