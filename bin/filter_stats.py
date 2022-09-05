#!/usr/bin/env python3
# @author Daniel Straub
# Takes two TSV count table from QIIME2
# and reports how much counts were filtered.

import pandas as pd
import sys

# argument check
if len(sys.argv) != 3:
    exit("Usage: count_table_max_reads.py <unfiltered_feature-table.tsv> <filtered_feature-table.tsv>")

# read tsv and skip first two rows
data_unfiltered = pd.read_csv(sys.argv[1], sep="\t", skiprows=None)  # count table
data_filtered = pd.read_csv(sys.argv[2], sep="\t", skiprows=[0])  # count table

# drop feature ids
df_unfiltered = data_unfiltered.drop(data_unfiltered.columns[0], axis=1)
df_filtered = data_filtered.drop(data_filtered.columns[0], axis=1)

# make sample count sums
sums_unfiltered = df_unfiltered.sum()
sums_filtered = df_filtered.sum()

# merge dataframes
out = sums_unfiltered.to_frame(name="unfiltered").join(sums_filtered.to_frame(name="filtered"))
out["lost"] = out["unfiltered"] - out["filtered"]
out["retained_percent"] = out["filtered"] / out["unfiltered"] * 100
out["lost_percent"] = 100 - out["retained_percent"]

# add column with sample names at beginning
out = out.rename_axis("sample").reset_index()

# rename columns
out = out.rename(columns={"unfiltered": "input_tax_filter", "filtered": "filtered_tax_filter"})

# write file
out.to_csv("count_table_filter_stats.tsv", sep="\t", index=False)
