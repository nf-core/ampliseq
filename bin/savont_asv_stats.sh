#!/bin/bash

awk '
    NR == 1 {
        # Extract sample names from header (skip first column)
        for (i = 2; i <= NF; i++) {
            samples[i-1] = $i
        }
        num_samples = NF - 1
        next
    }
    {
        # Accumulate counts for each sample
        for (i = 2; i <= NF; i++) {
            totals[samples[i-1]] += $i
        }
    }
    END {
        print "sample\tsavont_output"
        # Print in order of first appearance
        for (i = 1; i <= num_samples; i++) {
            print samples[i] "\t" totals[samples[i]]
        }
    }
' "$1"
