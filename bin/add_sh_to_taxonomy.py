#!/usr/bin/env python3
#@author Jeanette Tångrot

# Takes 

# Input: * ASV_tax.tsv
#        * "blastout.pb363.full0partial.tab"
#        * seq2sh.tsv
#        * SHs.tax
# Output: * ASV_tax_sh.tsv
# Usage: add_sh_to_taxonomy.py <seq2sh.tsv> <SHs.tax> <tax.tsv> <blastout.tab> <outfile>
#          <seq2sh.tsv>
#          <SHs.tax>
#          <tax.tsv>
#          <blastout.tab>
#          <outfile>

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
taxtable = pd.read_csv(sys.argv[3], sep='\t', header=0) #count table
# Old header line:
# ASV_ID  Domain  Kingdom Phylum  Class   Order   Family  Genus   confidence      sequence
# New header: 
# ASV_ID  Domain  Kingdom Phylum  Class   Order   Family  Genus  SH confidence      sequence
#DataFrame.insert(loc, column, value, allow_duplicates=False)
taxtable.insert(num_ranks+1,"SH","", allow_duplicates=False)
#print("taxtable:")
#print(taxtable)

# Go through vsearch matches and update taxonomy for those entries
fh = open( sys.argv[4], mode = 'r' )
prev_ASV = fh.readline().split()[0]
fh.seek(0)
matches = []
maxid = -1
maxlen = -1
for row in fh:
    #print(row)
    #[ASV, match, pid, alen, mism, gap_open, q_start, q_end, s_start, s_end, e_value, bit_score] = row.split()
    [ASV, match, pid, alen, therest] = row.split(maxsplit=4)
    pid = float(pid)
    alen = float(alen)
    
    if ASV != prev_ASV:
        SH = ""
        tax = ""
        conf = 0.0
        for m in matches:
#  >Nectria_magnispora|JF832665|SH1546528.08FU|refs|Fungi;Ascomycota;Sordariomycetes;Hypocreales;Nectriaceae;Nectria
            matchparts = m[0].split('|')
            #print( "Checking " +  matchparts[1] + "\n")
            try:
                new_SH = seq2sh.loc[ matchparts[1] ][1]
            except KeyError:
                print( "WARNING: " + matchparts[1] + " not in seq2SH list", file=sys.stderr )
                new_SH = ""
            if ( pd.isna( new_SH ) ) :
                print( "WARNING: no SH reported for " + matchparts[1], file=sys.stderr )
                new_SH = ""                
            #print( "   Checked " +  matchparts[1] + "; found **" + new_SH + "**\n")
            if SH != "" and new_SH != SH :
                SH = ""
                tax = ""
                break
            elif new_SH != "":
                SH = new_SH
                #tax = matchparts[4]
                #print( "SH: " + SH + "\n")
                try:
                    tax = list(shtax.loc[ SH ])
                except KeyError:
                    print( "WARNING: no taxonomy found for " + SH, file=sys.stderr )
                    tax = [""]*num_ranks
                #print( "tax: " + str(tax) + "\n")
                #conf = round( m[1]/100.0, 3 )
                conf = m[1]/100.0
                #print("conf " + str(conf) + " pid " + str(m[1]))
        #print( "SH: **" + SH + "**\n")
        if SH != "":
            #tax_list = tax.split(';') + [""]*(7 - len(tax.split(';'))) + [SH] + [conf]
            tax_list = tax[0:num_ranks] + [SH] + [conf]
            #print(tax_list)
            taxtable.loc[ taxtable['ASV_ID'] == prev_ASV, tax_entries] = tax_list
            #print("taxtable2:")
            #print(taxtable)
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

#print( "Last row;\n  " + row )
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
            #tax = matchparts[4]
            conf = m[1]/100.0
            #print("conf " + str(conf) + " pid " + str(m[1]))
    if SH != "":
        tax_list = tax[0:num_ranks] + [SH] + [conf]
        #print(tax_list)
        #tax_list = tax.split(';') + [""]*(7 - len(tax.split(';'))) + [SH] + [conf]
        taxtable.loc[ taxtable['ASV_ID'] == prev_ASV, tax_entries] = tax_list

   
# Write new taxtable, with SH and new taxonomy added if found
taxtable.to_csv(outfile, sep="\t", na_rep="", float_format="%.2f", index=False)

