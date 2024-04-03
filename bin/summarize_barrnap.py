#!/usr/bin/env python3
# @author Jeanette Tångrot
# Takes a list of files with barrnap predictions (rrna.arc.gff, rrna.bac.gff, etc)
# for ASV sequences, extracts evalues for each prediction and summarize the results
# in a new file "summary.gff". Assumes that the same program/barrnap version is
# used for all predictions. If there is more than one gene for a given domain,
# retains the lowest e-value (case of full rRNA operon sequences).

# import pandas as pd
import sys

# Initialize
method = dict()
evalues = dict()
orgs = set()

# Go through each file and store evalues for all predictions for each query sequence
for file in sys.argv[1:]:
    org = file.removeprefix("rrna.")
    org = org.replace(".gff", "_eval")
    orgs.add(org)
    fh = open(file, mode="r")
    for row in fh:
        if row.startswith("#"):
            continue
        rowparts = row.split()
        asv = rowparts[0]
        method[asv] = rowparts[1]
        if asv not in evalues:
            evalues[asv] = dict()
        if (org not in evalues[asv]) or (float(evalues[asv][org]) > float(rowparts[5])):
           evalues[asv][org] = rowparts[5]
    fh.close()

# Write results
fh = open("summary.tsv", mode="w")
orglist = sorted(orgs)
header = sorted(orgs)
header.insert(0, "ASV_ID")
header.append("eval_method")
fh.write("\t".join(header) + "\n")
for asv, meth in method.items():
    row = [asv]
    for org in orglist:
        if org in evalues[asv]:
            row.append(evalues[asv][org])
        else:
            row.append("NA")
    row.append(meth)
    fh.write("\t".join(row) + "\n")
fh.close()
