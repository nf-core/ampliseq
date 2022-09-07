#!/usr/bin/env Rscript

# sbdiexportreannotate.R
#
# A script that collates taxonomy data from Ampliseq to produce an updated
# annotation tsv file as close to ready for submission to the Swedish
# Biodiversity Data Infrastructure (SBDI) as possible.
#
# The script expects the following arguments: dbversion, ASV_tax_species.tsv
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(tidyverse))

# Get dbversion and taxonomy file from the command line
args            <- commandArgs(trailingOnly=TRUE)
dbversion       <- args[1]
taxfile         <- args[2]

taxonomy <- read.delim(taxfile, sep = '\t', stringsAsFactors = FALSE)

taxonomy %>%
    mutate(SH = if("SH" %in% colnames(.)) SH else '') %>%
    relocate(SH, .after = Species) %>%
    rename_with(tolower, Domain:Species) %>%
    rename(
        asv_id_alias = ASV_ID,
        asv_sequence = sequence,
        specificEpithet = species,
        otu = SH,
        annotation_confidence = confidence
    ) %>%
    mutate(
        infraspecificEpithet = '',
        domain = str_remove(domain, 'Reversed:_'),
        scientificName = case_when(
            !(is.na(specificEpithet) | specificEpithet == '') ~ sprintf("%s %s", genus, specificEpithet),
            !(is.na(genus)   | genus == '')                   ~ sprintf("%s", genus),
            !(is.na(family)  | family == '')                  ~ sprintf("%s", family),
            !(is.na(order)   | order == '')                   ~ sprintf("%s", order),
            !(is.na(class)   | class == '')                   ~ sprintf("%s", class),
            !(is.na(phylum)  | phylum == '')                  ~ sprintf("%s", phylum),
            !(is.na(kingdom) | kingdom == '')                 ~ sprintf("%s", kingdom),
            TRUE                                              ~ 'Unassigned'
        ),
        taxonRank = case_when(
            !(is.na(specificEpithet) | specificEpithet == '') ~ 'species',
            !(is.na(genus)   | genus == '')                   ~ 'genus',
            !(is.na(family)  | family == '')                  ~ 'family',
            !(is.na(order)   | order == '')                   ~ 'order',
            !(is.na(class)   | class == '')                   ~ 'class',
            !(is.na(phylum)  | phylum == '')                  ~ 'phylum',
            !(is.na(kingdom) | kingdom == '')                 ~ 'kingdom',
            TRUE                                              ~ 'kingdom'
        ),
        date_identified = as.character(lubridate::today()),
        reference_db = dbversion,
        annotation_algorithm = case_when(
            (!(is.na(otu) | otu == '')) ~ 'Ampliseq:addsh',
            TRUE                        ~ 'DADA2:assignTaxonomy:addSpecies'
        ),
        identification_references = 'https://docs.biodiversitydata.se/analyse-data/molecular-tools/#taxonomy-annotation',
        taxon_remarks = '',
        kingdom = ifelse(is.na(kingdom), 'Unassigned', kingdom)
    ) %>%
    relocate(asv_sequence, .after = asv_id_alias) %>%
    relocate(infraspecificEpithet:identification_references, .after = specificEpithet) %>%
    relocate(otu, .after = infraspecificEpithet) %>%
    select(-domain) %>%
    write_tsv("annotation.tsv", na = '')
