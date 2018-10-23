#!/usr/bin/env python3
#@author Daniel Straub
# Takes two CSV files from QIIME2 demux output and a quality threshold 
# to generate a tuple of index locations that resemble the cutoff value used for DADA2 in QIIME2.

import pandas as pd
import sys 

#argument check
if len(sys.argv) != 4:
    exit("Usage: dada_trunc_parameter.py <forward-seven-number-summaries.csv> <reverse-seven-number-summaries.csv> <int>")

#parameters
data_fw = pd.read_csv(sys.argv[1]) #quality values forward reads
data_rv = pd.read_csv(sys.argv[2]) #quality values reverse reads
qmin = float(sys.argv[3]) #threshold

#select row with median values (file row 6, starting with "50%") and drop first row
median_fw = data_fw.iloc[4][1:].values.tolist()
median_rv = data_rv.iloc[4][1:].values.tolist()

#iterate through values and find first value that falls below threshold
def function(median):
    trunc = len(median)
    for value in median:
        if value < qmin:
            trunc = median.index(value)
            break
    return trunc
trunc_fw = function(median_fw)
trunc_rv = function(median_rv)

#print values
print(trunc_fw, trunc_rv, sep=',', end='')