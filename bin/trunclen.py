#!/usr/bin/env python3
# @author Daniel Straub
# Takes two CSV files from QIIME2 demux output, a quality threshold and a cutoff for the retained read fraction
# to generate a tuple of index locations that resemble the cutoff value used for DADA2 in QIIME2.

import pandas as pd
import sys

# argument check
if len(sys.argv) != 4:
    exit("Usage: dada_trunc_parameter.py <*_qual_stats.tsv> <int[0-40]> <float[0-1]>")

# parameters
data = pd.read_csv(sys.argv[1], delimiter="\t")  # quality values forward reads
qmin = float(sys.argv[2])  # quality threshold
rmin = float(sys.argv[3])  # read count threshold (fraction)

# select row with median values (file row 6, starting with "50%") and drop first row
median = data.iloc[1][1:].values.tolist()

# select row with count numbers (file row name "count")
reads = data.iloc[0][1:].values.tolist()
# extract maximum read count
fraction_reads = int(max(reads) * rmin)


# iterate through values and find first value that falls below threshold
def function(values, cutoff):
    trunc = len(values)
    for value in values:
        if value < cutoff:
            trunc = values.index(value)
            break
    return trunc


# find quality threshold
trunc_median = function(median, qmin)

# find read threshold
trunc_reads = function(reads, fraction_reads)

# final threshold
trunc = min(trunc_median, trunc_reads)

# print values
print(trunc, end="")
