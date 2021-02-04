#!/usr/bin/env python3
#
# Takes two TSV files as input; stats from dada2 filtering and stats
# from dada2 denoising, and reports these numbers together with
# fractions of input sequences that are filtered/non-chimeric.
# Results are written to file "dada_stats.tsv"
#
# Jeanette Tångrot

#-- Import libraries, do initializations  --#
import pandas as pd
import sys 

file_out = "dada_stats.tsv"

#-- Check arguments --#
if len( sys.argv ) != 3:
    exit( "Usage: make_dada2_stats_pacbio.py <filter_stats.tsv> <denoise_stats.tsv>" )

#-- Reas TSVs --#
filt = pd.read_csv( sys.argv[1], sep = '\t', usecols = ['sample.id', 'file', 'reads.in', 'reads.out'] )
denoise = pd.read_csv( sys.argv[2], sep = '\t' )

#-- Count number of input sequences --#
num_input = filt[ 'reads.in' ]

#-- Create results table --#
res = filt.join( denoise.set_index( 'file' ), on = 'file' )
res.pop( 'file' )
res['perc_filt'] = res['reads.out'] / num_input * 100
res['perc_denoise'] = res['denoised'] / num_input * 100
res['perc_nonchim'] = res['nonchim'] / num_input * 100
res = res[ ['sample.id', 'reads.in', 'reads.out', 'perc_filt', 'denoised', 'perc_denoise', 'nonchim', 'perc_nonchim'] ]
res.to_csv( file_out, sep = '\t', header = [ 'sample-id', 'input', 'filtered', 'percentage of input passed filter', 'denoised', 'percentage of input denoised', 'non-chimeric', 'percentage of input non-chimeric' ], index=False )

