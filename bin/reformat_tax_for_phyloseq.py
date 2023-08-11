#!/usr/bin/env python3

import pandas as pd
import sys

tax_file = sys.argv[1]
out_file = sys.argv[2]

# Import tsv file
tax_df = pd.read_csv(tax_file, sep="\t")

# The second column should hold the taxonomy information
tax_col = tax_df.columns[1]

# Split the values in the tax column
split_tax = tax_df[tax_col].str.split(";", expand=True)

# Assign names to the new columns with an auto incrementing integer
new_col_names = [f"{tax_col}_{i+1}" for i in range(split_tax.shape[1])]
split_tax.columns = new_col_names

# Strip whitespace from the tax names
split_tax = split_tax.applymap(lambda x: x.strip() if isinstance(x, str) else x)

# Drop the original tax column
tax_df = tax_df.drop(columns=[tax_col])

# Add the new tax columns to the df
result = pd.concat([tax_df, split_tax], axis=1)

# Create new tsv file
result.to_csv(out_file, sep="\t", index=False)
