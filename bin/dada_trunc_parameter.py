#!/usr/bin/env python3
#@author Daniel Straub
# Takes two CSV files from QIIME2 demux output, a quality threshold and a cutoff for the retained read fraction
# to generate a tuple of index locations that resemble the cutoff value used for DADA2 in QIIME2.

import pandas as pd
import sys 

#argument check
if len(sys.argv) != 5:
    exit("Usage: dada_trunc_parameter.py <forward-seven-number-summaries.csv> <reverse-seven-number-summaries.csv> <int[0-40]> <float[0-1]>")

#parameters
data_fw = pd.read_csv(sys.argv[1]) #quality values forward reads
data_rv = pd.read_csv(sys.argv[2]) #quality values reverse reads
qmin = float(sys.argv[3]) #quality threshold
rmin = float(sys.argv[4]) #read count threshold (fraction)

#select row with median values (file row 6, starting with "50%") and drop first row
median_fw = data_fw.iloc[4][1:].values.tolist()
median_rv = data_rv.iloc[4][1:].values.tolist()

#select row with count numbers (file row name "count")
reads_fw = data_fw.iloc[0][1:].values.tolist()
reads_rv = data_rv.iloc[0][1:].values.tolist()
#extract maximum read count
fraction_reads = int(max(reads_fw)*rmin)

#iterate through values and find first value that falls below threshold
def function(values, cutoff):
    trunc = len(values)
    for value in values:
        if value < cutoff:
            trunc = values.index(value)
            break
    return trunc

#find quality threshold
trunc_median_fw = function(median_fw, qmin)
trunc_median_rv = function(median_rv, qmin)

#find read threshold
trunc_reads_fw = function(reads_fw, fraction_reads)
trunc_reads_rv = function(reads_rv, fraction_reads)

#final threshold
trunc_fw = min(trunc_median_fw,trunc_reads_fw)
trunc_rv = min(trunc_median_rv,trunc_reads_rv)

#print values
print(trunc_fw, trunc_rv, sep=',', end='')