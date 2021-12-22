#!/usr/bin/env Rscript

# sbdiexportreannotate.R
#
# A script that collates taxonomy data from Ampliseq to produce an Excel file as close
# to ready for submission to the Swedish Biodiversity Data Infrastructure
# (SBDI) as possible.
#
# The script expects the following positional arguments:
#   1. dbversion
#   2. Input table, e.g. ASV_tax_species.tsv.
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(tidyverse))

# Get the library layout and primers from the command line
args            <- commandArgs(trailingOnly=TRUE)
dbversion       <- args[1]
taxfile         <- args[2]

taxonomy <- read.delim(taxfile, sep = '\t', stringsAsFactors = FALSE) %>%
    # Pick ranks by number, to circumvent that PR2 has different ranks than the rest
    rename(Domain = 2, Kingdom = 3, Phylum = 4, Class = 5, Order = 6, Family = 7, Genus = 8, Species = 9)

taxonomy %>%
    # I'm keeping this if the ugly renaming above proves unnecessary one day
    rename_with(tolower, Domain:Species) %>%
    rename(
        asv_id_alias = ASV_ID,
        asv_sequence = sequence,
        specificEpithet = species,
        annotation_confidence = confidence
    ) %>%
    mutate(
        infraspecificEpithet = '',
        otu = '',
        domain = str_remove(domain, 'Reversed:_'),
        scientificName = case_when(
            !is.na(specificEpithet) ~ sprintf("%s %s", genus, specificEpithet),
            !is.na(genus)           ~ sprintf("%s", genus),
            !is.na(family)          ~ sprintf("%s", family),
            !is.na(order)           ~ sprintf("%s", order),
            !is.na(class)           ~ sprintf("%s", class),
            !is.na(phylum)          ~ sprintf("%s", phylum),
            !is.na(kingdom)          ~ sprintf("%s", kingdom),
            TRUE                    ~ 'Unassigned'
        ),
        taxonRank = case_when(
            !is.na(specificEpithet) ~ 'species',
            !is.na(genus)           ~ 'genus',
            !is.na(family)          ~ 'family',
            !is.na(order)           ~ 'order',
            !is.na(class)           ~ 'class',
            !is.na(phylum)          ~ 'phylum',
            !is.na(kingdom)         ~ 'kingdom',
            TRUE                    ~ 'kingdom'
        ),
        date_identified = as.character(lubridate::today()),
        reference_db = dbversion,
        annotation_algorithm = 'DADA2:assignTaxonomy:addSpecies',
        identification_references = 'https://docs.biodiversitydata.se/analyse-data/molecular-tools/#taxonomy-annotation',
        taxon_remarks = '',
        kingdom = ifelse(is.na(kingdom), 'Unassigned', kingdom)
    ) %>%
    relocate(asv_sequence, .after = asv_id_alias) %>%
    relocate(infraspecificEpithet:identification_references, .after = specificEpithet) %>%
    select(-domain) %>%
    write_tsv("annotation.tsv", na = '')
