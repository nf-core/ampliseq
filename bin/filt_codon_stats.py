#!/usr/bin/env python3

import pandas as pd
import sys

# argument check
if len(sys.argv) != 2:
    exit("Usage: count_table_max_reads.py <ASVcounts.tsv>")

# read tsv
data = pd.read_csv(sys.argv[1], sep="\t", skiprows=None)  # count table

# drop feature ids
df = data.drop(data.columns[0], axis=1)

# sum each column
sums = df.sum(axis=0)

# add column with sample names at beginning
out = sums.rename_axis("sample").reset_index()

# rename columns
out = out.rename(columns={0: "encoding"})

# write file
out.to_csv("codon.filtered.stats.tsv", sep="\t", index=False)
