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
  'event_id_alias' = colnames(asvs)[2:(n_samples+1)],
  'materialSampleID' = character(n_samples),
  'eventDate' = character(n_samples),
  'samplingProtocol' = character(n_samples),
  'locationID' = character(n_samples),
  'decimalLatitude' = character(n_samples),
  'decimalLongitude' = character(n_samples),
  'geodeticDatum' = character(n_samples),
  'coordinateUncertaintyInMeters' = character(n_samples),
  'recordedBy' = character(n_samples),
  'country' = character(n_samples),
  'municipality' = character(n_samples),
  'verbatimLocality' = character(n_samples),
  'minimumElevationInMeters' = character(n_samples),
  'maximumElevationInMeters' = character(n_samples),
  'minimumDepthInMeters' = character(n_samples),
  'maximumDepthInMeters' = character(n_samples)
) %>%
  inner_join(
    asvs %>% pivot_longer(2:ncol(.), names_to = 'event_id_alias', values_to = 'count') %>%
      group_by(event_id_alias) %>% summarise(sampleSizeValue = sum(count), .groups = 'drop'),
    by = 'event_id_alias'
  ) %>%
  select(1:4, sampleSizeValue, 5:17) %>%
  write_tsv("event.tsv", na = '')

# mixs
data.frame(
  'event_id_alias' = colnames(asvs)[2:(n_samples+1)],
  'sop' = character(n_samples),
  'target_gene' = character(n_samples),
  'target_subfragment' = character(n_samples),
  'lib_layout' = character(n_samples),
  'pcr_primer_name_forward' = character(n_samples),
  'pcr_primer_name_reverse' = character(n_samples),
  'pcr_primer_forward' = character(n_samples),
  'pcr_primer_reverse' = character(n_samples),
  'env_broad_scale' = character(n_samples),
  'env_local_scale' = character(n_samples),
  'env_medium' = character(n_samples)
) %>%
  write_tsv("mixs.tsv", na = '')

# emof
data.frame(
  'event_id_alias' = character(n_samples),
  'measurementType' = character(n_samples),
  'measurementTypeID' = character(n_samples),
  'measurementValue' = character(n_samples),
  'measurementValueID' = character(n_samples),
  'measurementUnit' = character(n_samples),
  'measurementUnitID' = character(n_samples),
  'measurementAccuracy' = character(n_samples),
  'measurementDeterminedDate' = character(n_samples),
  'measurementDeterminedBy' = character(n_samples),
  'measurementMethod' = character(n_samples),
  'measurementRemarks' = character(n_samples)
) %>%
  write_tsv("emof.tsv", na = '')

# asv-table
asvtax   <- asvs %>% 
  inner_join(taxonomy, by = 'ASV_ID') %>%
  rename(asv_id_alias = ASV_ID, DNA_sequence = sequence) %>%
  rename_with(tolower, Domain:Species) %>%
  rename(specificEpithet = species) %>%
  select(-confidence) %>%
  mutate(associatedSequences = '', domain = str_remove(domain, 'Reversed:_')) %>%
  write_tsv("asv-table.tsv", na = '')
