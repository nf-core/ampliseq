#!/usr/bin/env python3

import argparse
import gzip
import pandas as pd
import sys

usage = """This program filters ASVs that aren't centroids after post-clustering."""

parser = argparse.ArgumentParser(description=usage)

parser.add_argument(
    "-t",
    "--count-table",
    dest="count",
    type=argparse.FileType("r"),
    help="Count table file of the ASVs from the AmpliSeq pipeline",
    required=True,
)

parser.add_argument(
    "-p",
    "--prefix",
    dest="prefix",
    type=str,
    help="Prefix of the output files",
    required=True,
)

parser.add_argument(
    "-c",
    "--cluster-fastas",
    dest="cluster_fastas",
    type=argparse.FileType('r'),
    default=sys.stdin,
    help="Space separated list of fasta files of the clusters. First read of the cluster should be the centroid of that cluster.",
    required=True,
)

count = parser.parse_args().count
prefix = parser.parse_args().prefix

# This dictionary will store the centroid ASVs as keys, and the values will be the ASVs clustered to that centroid
cluster_dict = {}

# Loop though list of cluster fasta files to populate cluster_dict and to create centroid fasta file
cluster_fastas = parser.parse_args().cluster_fastas.read().rstrip().split(" ")
for cluster_fasta in cluster_fastas:
    read_num = 0

    # Loop through each line of current fasta file and open output fasta file in append mode
    with gzip.open(cluster_fasta, "rt") as in_fasta, open(prefix + "_filtered.fna", "a") as out_fasta:
        for line in in_fasta:
            line = line.rstrip("\n")

            # If the line is not a sequence
            if line.startswith(">"):
                read_num += 1
                asv_name = line[1:]

                # If the read is the centroid
                if read_num == 1:
                    centroid_name = asv_name
                    cluster_dict[centroid_name] = []
                    out_fasta.write(f"{line}\n")

                # If the read is not the centroid
                else:
                    cluster_dict[centroid_name].append(asv_name)

            # If the line is a sequence
            else:
                # If the read is the centroid
                if read_num == 1:
                    out_fasta.write(f"{line}\n")

# This dictionary will store the samples as keys, and the values will be the number of ASVs in that sample
sam_asv_counts = {}

# This count_df will have ASVs as the index, and samples as the header
count_df = pd.read_table(count, delimiter="\t", index_col=0, header=0)

# Get the number of ASVs per sample before clustering
for sample in count_df.columns:
    sam_asv_counts[sample] = (count_df[sample] != 0).sum()
stats_df = pd.DataFrame(list(sam_asv_counts.items()), columns=["sample", "ASVs_before_clustering"])

# Loop through centroids
for centroid in cluster_dict.keys():
    # If the current centroid has ASVs clustered to it
    if cluster_dict[centroid] != []:
        # Get a list of all ASVs in the cluster (including the centroid)
        cluster_list = cluster_dict[centroid].copy()
        cluster_list.append(centroid)

        # Sum all rows in the cluster
        summed_row = count_df.loc[cluster_list].sum()

        # Assign summed row to centroid row and drop non-centroid rows
        count_df.loc[centroid] = summed_row
        count_df.drop(cluster_dict[centroid], inplace=True)

# Get the number of ASVs per sample after clustering
for sample in count_df.columns:
    sam_asv_counts[sample] = (count_df[sample] != 0).sum()
stats_df["ASVs_after_clustering"] = list(sam_asv_counts.values())

# Output filtered count tsv and stats tsv
count_df.to_csv(prefix + "_filtered.table.tsv", sep="\t")
stats_df.to_csv(prefix + "_filtered.stats.tsv", sep="\t", index=False)
