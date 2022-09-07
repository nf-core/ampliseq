#!/usr/bin/env python3
# @author Daniel Straub
# Takes one TSV count table from QIIME2
# and reports the maximum or minimum counts of all samples.

import pandas as pd
import sys

# argument check
if len(sys.argv) != 3 or sys.argv[2] not in ["maximum", "minimum"]:
    exit("Usage: count_table_max_reads.py <feature-table.tsv> <maximum/minimum>")

# read tsv and skip first two rows
data = pd.read_csv(sys.argv[1], sep="\t", skiprows=[0, 1], header=None)  # count table

# drop feature ids
df = data.drop(data.columns[0], axis=1)

# make sums
sums = df.sum()

# determine maximum or minimum
if sys.argv[2] == "maximum":
    out = int(sums.max())
elif sys.argv[2] == "minimum":
    out = int(sums.min())

# print value
print(out, end="")
