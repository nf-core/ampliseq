#!/usr/bin/env Rscript

# sbdiexport.R
#
# A script that collates data from Ampliseq to produce an Excel file as close 
# to ready for submission to the Swedish Biodiversity Data Infrastructure
# (SBDI) as possible.
#
# The script expects the following input files to be present in the directory:
# ASV_table.tsv and ASV_tax_species.tsv.
#
# It also expects the following arguments: paired|single fwdprimer revprimer
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(tidyverse))

# Basename for outfiles
BASENAME = "SBDI-metabarcoding-template"

# Missing params that we might ask for, otherwise provide default values
sop                     <- "STANDARD OPERATING PROCEDURE"
target_gene             <- "TARGET GENE"
target_subfragment      <- ""
pcr_primer_name_forward <- "FWD PRIMER NAME"
pcr_primer_name_reverse <- "REV PRIMER NAME"
env_broad_scale         <- "SEE MIXS SPEC"
env_local_scale         <- "SEE MIXS SPEC"
env_medium              <- "SEE MIXS SPEC"
event_id_alias          <- ""
measurementType         <- ""
measurementTypeID       <- ""
measurementValue        <- ""
measurementValueID      <- ""
measurementUnit         <- ""
measurementUnitID       <- ""
measurementAccuracy     <- ""
measurementDeterminedDate <- ""
measurementDeterminedBy <- ""
measurementMethod       <- ""
measurementRemarks      <- ""

# Get the library layout and primers from the command line
args = commandArgs(trailingOnly=TRUE)
lib_layout      <- args[1]
fwd_primer_seq  <- args[2]
rev_primer_seq  <- args[3]

write(sprintf("DEBUG: args: %s", paste(args, collapse = ", ")), stderr())

asvs <- read.delim("ASV_table.tsv", sep = '\t')
n_samples <- length(colnames(asvs)) - 1

taxonomy <- read.delim("ASV_tax_species.tsv", sep = '\t')

# Write the event tab
data.frame(
  'event_id_alias' = colnames(asvs)[2:(n_samples+1)], 'materialSampleID' = character(n_samples),
  'eventDate' = character(n_samples), 'samplingProtocol' = character(n_samples),
  'sampleSizeValue' = character(n_samples), 'locationID' = character(n_samples),
  'decimalLatitude' = character(n_samples), 'decimalLongitude' = character(n_samples),
  'geodeticDatum' = character(n_samples), 'coordinateUncertaintyInMeters' = character(n_samples),
  'recordedBy' = character(n_samples), 'country' = character(n_samples),
  'municipality' = character(n_samples), 'verbatimLocality' = character(n_samples),
  'minimumElevationInMeters' = character(n_samples), 'maximumElevationInMeters' = character(n_samples),
  'minimumDepthInMeters' = character(n_samples), 'maximumDepthInMeters' = character(n_samples)
) %>%
  write_tsv(sprintf("%s.event.tsv", BASENAME))

# mixs
data.frame(
  'event_id_alias' = colnames(asvs)[2:(n_samples+1)],
  'sop' = rep(sop, n_samples),
  'target_gene' = rep(target_gene, n_samples),
  'target_subfragment' = rep(target_subfragment, n_samples),
  'lib_layout' = rep(lib_layout, n_samples), 
  'pcr_primer_seq_name_forward' = rep(pcr_primer_name_forward, n_samples),
  'pcr_primer_seq_name_reverse' = rep(pcr_primer_name_reverse, n_samples),
  'pcr_primer_seq_forward' = rep(fwd_primer_seq, n_samples),
  'pcr_primer_seq_reverse' = rep(rev_primer_seq, n_samples), 
  'env_broad_scale' = rep(env_broad_scale, n_samples),
  'env_local_scale' = rep(env_local_scale, n_samples),
  'env_medium' = rep(env_medium, n_samples)
) %>%
  write_tsv(sprintf("%s.mixs.tsv", BASENAME))

# emof
data.frame(
  'event_id_alias' = colnames(asvs)[2:(n_samples+1)],
  'measurementType' = rep(measurementType, n_samples),
  'measurementTypeID' = rep(measurementTypeID, n_samples),
  'measurementValue' = rep(measurementValue, n_samples),
  'measurementValueID' = rep(measurementValueID, n_samples),
  'measurementUnit' = rep(measurementUnit, n_samples),
  'measurementUnitID' = rep(measurementUnitID, n_samples),
  'measurementAccuracy' = rep(measurementAccuracy, n_samples),
  'measurementDeterminedDate' = rep(measurementDeterminedDate, n_samples),
  'measurementDeterminedBy' = rep(measurementDeterminedBy, n_samples),
  'measurementMethod' = rep(measurementMethod, n_samples),
  'measurementRemarks' = rep(measurementRemarks, n_samples)
) %>%
  write_tsv(sprintf("%s.emof.tsv", BASENAME))

# asv-table
asvtax   <- asvs %>% 
  inner_join(taxonomy, by = 'ASV_ID') %>%
  rename(asv_id_alias = ASV_ID, DNA_sequence = sequence) %>%
  rename_with(tolower, Domain:Species) %>%
  write_tsv(sprintf("%s.asv-table.tsv", BASENAME))
