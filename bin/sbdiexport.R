#!/usr/bin/env Rscript

# sbdiexport.R
#
# A script that collates data from Ampliseq to produce four tsv files as close
# as possible to ready for submission to the Swedish Biodiversity Data
# Infrastructure (SBDI) as possible.
#
# The script expects the following arguments: paired|single fwdprimer revprimer asvtable taxonomytable [metadata]
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(tidyverse))

EVENT_COLS <- c(
    'materialSampleID', 'eventDate', 'samplingProtocol', 'locationID', 'decimalLatitude',
    'decimalLongitude', 'geodeticDatum', 'coordinateUncertaintyInMeters', 'recordedBy', 'country',
    'municipality', 'verbatimLocality', 'minimumElevationInMeters', 'maximumElevationInMeters',
    'minimumDepthInMeters', 'maximumDepthInMeters'
)
MIXS_COLS  <- c(
    'sop', 'target_gene', 'target_subfragment', 'pcr_primer_name_forward',
    'pcr_primer_name_reverse', 'env_broad_scale', 'env_local_scale', 'env_medium'
)
EMOF_COLS  <- c(
    'measurementType', 'measurementTypeID', 'measurementValue', 'measurementValueID',
    'measurementUnit', 'measurementUnitID', 'measurementAccuracy', 'measurementDeterminedDate',
    'measurementDeterminedBy', 'measurementMethod', 'measurementRemarks'
)

# Get the library layout and primers from the command line
args            <- commandArgs(trailingOnly=TRUE)
lib_layout      <- args[1]
fwd_primer_seq  <- args[2]
rev_primer_seq  <- args[3]
asvtable        <- args[4]
taxtable        <- args[5]
metadata        <- args[6]

asvs <- read.delim(asvtable, sep = '\t', stringsAsFactors = FALSE)
n_samples <- length(colnames(asvs)) - 1

taxonomy <- read.delim(taxtable, sep = '\t', stringsAsFactors = FALSE)

# Read the metadata table if provided, otherwise create one
if ( ! is.na(metadata) ) {
    meta <- read.delim(metadata, sep = '\t', stringsAsFactors = FALSE)
} else {
    meta <- data.frame(ID = colnames(asvs)[2:(n_samples+1)])
}

# Make sure it's congruent with the asv table
meta <- data.frame(
    'event_id_alias' = colnames(asvs)[2:(n_samples+1)]
) %>%
    left_join(meta, by = c('event_id_alias' = 'ID')) %>%
    distinct(event_id_alias, .keep_all = TRUE) %>%
    arrange(event_id_alias)

# Write the event tab
event <- data.frame(
    'event_id_alias' = colnames(asvs)[2:(n_samples+1)]
) %>%
    arrange(event_id_alias)
for ( c in EVENT_COLS ) {
    if ( c %in% colnames(meta) ) {
        event[[c]] <- meta[[c]]
    } else {
        event[[c]] <- character(n_samples)
    }
}
event %>%
    inner_join(
        asvs %>% pivot_longer(2:ncol(.), names_to = 'event_id_alias', values_to = 'count') %>%
        group_by(event_id_alias) %>% summarise(sampleSizeValue = sum(count), .groups = 'drop'),
        by = 'event_id_alias'
    ) %>%
    select(1:4, sampleSizeValue, 5:ncol(.)) %>%
    write_tsv("event.tsv", na = '')

# mixs
mixs <- data.frame(
    'event_id_alias' = colnames(asvs)[2:(n_samples+1)],
    'lib_layout' = rep(lib_layout, n_samples),
    'pcr_primer_forward' = rep(fwd_primer_seq, n_samples),
    'pcr_primer_reverse' = rep(rev_primer_seq, n_samples)
) %>%
    arrange(event_id_alias)
for ( c in MIXS_COLS ) {
    if ( c %in% colnames(meta) ) {
        mixs[[c]] <- meta[[c]]
    } else {
        mixs[[c]] <- character(n_samples)
    }
}
mixs %>%
    relocate(lib_layout, .after = target_subfragment) %>%
    relocate(pcr_primer_forward, pcr_primer_reverse, .after = pcr_primer_name_reverse) %>%
    write_tsv("mixs.tsv", na = '')

# emof
emof <- data.frame(
    'event_id_alias' = colnames(asvs)[2:(n_samples+1)]
) %>%
    arrange(event_id_alias)
for ( c in EMOF_COLS ) {
    if ( c %in% colnames(meta) ) {
        emof[[c]] <- meta[[c]]
    } else {
        emof[[c]] <- character(n_samples)
    }
}
emof %>% write_tsv("emof.tsv", na = '')

# asv-table
asvtax <- asvs %>%
    inner_join(taxonomy, by = 'ASV_ID') %>%
    rename_with(tolower, Domain:Species) %>%
    rename(
        specificEpithet = species,
        asv_id_alias = ASV_ID,
        DNA_sequence = sequence
    ) %>%
    mutate(
        domain = str_remove(domain, 'Reversed:_'),
        associatedSequences = '',
        infraspecificEpithet = '',
        otu = '',
        kingdom = ifelse(is.na(kingdom), 'Unassigned', kingdom)
    ) %>%
    relocate(DNA_sequence:associatedSequences, .before = domain) %>%
    select(-confidence, -domain) %>%
    write_tsv("asv-table.tsv", na = '')
