#!/usr/bin/env Rscript

# sbdiexport.R
#
# A script that collates data from Ampliseq to produce four tsv files as close
# to ready for submission to the Swedish Biodiversity Data
# Infrastructure (SBDI) as possible.
#
# The script expects the following arguments: paired|single fwdprimer revprimer asvtable taxonomytable [metadata]
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(tidyverse))

EVENT_COLS <- c(
    'datasetID', 'institutionCode', 'institutionID', 'collectionCode', 'materialSampleID',
    'associatedSequences', 'fieldNumber', 'catalogNumber', 'references', 'eventDate', 'samplingProtocol',
    'locationID', 'decimalLatitude', 'decimalLongitude', 'geodeticDatum', 'coordinateUncertaintyInMeters',
    'recordedBy', 'country', 'municipality', 'verbatimLocality', 'minimumElevationInMeters',
    'maximumElevationInMeters', 'minimumDepthInMeters', 'maximumDepthInMeters'
)
DNA_COLS  <- c(
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

# Read taxonomy table and make sure all expected columns are there
taxonomy <- read.delim(taxtable, sep = '\t', stringsAsFactors = FALSE) %>%
    mutate(Domain = if("Domain" %in% colnames(.)) Domain else '') %>%
    mutate(Kingdom = if("Kingdom" %in% colnames(.)) Kingdom else if ("Supergroup" %in% colnames(.)) Supergroup else '') %>%
    mutate(Phylum = if("Phylum" %in% colnames(.)) Phylum else if ("Division" %in% colnames(.)) Division else '') %>%
    mutate(Class = if("Class" %in% colnames(.)) Class else '') %>%
    mutate(Order = if("Order" %in% colnames(.)) Order else '') %>%
    mutate(Family = if("Family" %in% colnames(.)) Family else '') %>%
    mutate(Genus = if("Genus" %in% colnames(.)) Genus else '') %>%
    mutate(Species = if("Species" %in% colnames(.)) Species else '') %>%
    mutate(Species_exact = if("Species_exact" %in% colnames(.)) Species_exact else '') %>%
    mutate(otu = if("SH" %in% colnames(.)) SH else if ("BOLD_bin" %in% colnames(.)) BOLD_bin else '') %>%
    relocate(Domain, .after = sequence) %>%
    relocate(Kingdom, .after = Domain) %>%
    relocate(Phylum, .after = Kingdom) %>%
    relocate(Class, .after = Phylum) %>%
    relocate(Order, .after = Class) %>%
    relocate(Family, .after = Order) %>%
    relocate(Genus, .after = Family) %>%
    relocate(Species, .after = Genus) %>%
    relocate(Species_exact, .after = Species) %>%
    relocate(otu, .after = Species_exact)


# Read the metadata table if provided, otherwise create one
if ( ! is.na(metadata) ) {
    meta <- read.delim(metadata, sep = '\t', stringsAsFactors = FALSE)
} else {
    meta <- data.frame(ID = colnames(asvs)[2:(n_samples+1)])
}

# Make sure it's congruent with the asv table
meta <- data.frame(
    'eventID' = colnames(asvs)[2:(n_samples+1)]
) %>%
    left_join(meta, by = c('eventID' = 'ID')) %>%
    distinct(eventID, .keep_all = TRUE) %>%
    arrange(eventID)

# Write the event tab
event <- data.frame(
    'eventID' = colnames(asvs)[2:(n_samples+1)]
) %>%
    arrange(eventID)
for ( c in EVENT_COLS ) {
    if ( c %in% colnames(meta) ) {
        event[[c]] <- meta[[c]]
    } else {
        event[[c]] <- character(n_samples)
    }
}

# Add links to ENA
event$'materialSampleID' <- sub("^ERS", "https://www.ebi.ac.uk/ena/browser/view/ERS", event$'materialSampleID')
event$'materialSampleID' <- sub("^SAMEA", "https://www.ebi.ac.uk/ena/browser/view/SAMEA", event$'materialSampleID')
event$'associatedSequences' <- sub("^ERR", "https://www.ebi.ac.uk/ena/browser/view/ERR", event$'associatedSequences')

event %>%
    write_tsv("event.tsv", na = '')

# dna (previously mixs)
dna <- data.frame(
    'eventID' = colnames(asvs)[2:(n_samples+1)],
    'lib_layout' = rep(lib_layout, n_samples),
    'pcr_primer_forward' = rep(fwd_primer_seq, n_samples),
    'pcr_primer_reverse' = rep(rev_primer_seq, n_samples)
) %>%
    arrange(eventID)
for( c in DNA_COLS ) {
    if ( c %in% colnames(meta) ) {
        dna[[c]] <- meta[[c]]
    } else if ( c %in% c("sop") ) {
        dna[[c]] <- rep('https://nf-co.re/ampliseq',n_samples)
    } else {
        dna[[c]] <- character(n_samples)
    }
}
dna %>%
    relocate(lib_layout, .after = target_subfragment) %>%
    relocate(pcr_primer_forward, pcr_primer_reverse, .after = pcr_primer_name_reverse) %>%
    write_tsv("dna.tsv", na = '')

# emof
emof <- data.frame(
    'eventID' = colnames(asvs)[2:(n_samples+1)]
) %>%
    arrange(eventID)
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
    mutate(across(domain:species, ~str_replace_all(.,' ','_'))) %>%
    rename(
        specificEpithet = species,
        asv_id_alias = ASV_ID,
        DNA_sequence = sequence
    ) %>%
    mutate(
        domain = str_remove(domain, 'Reversed:_'),
        associatedSequences = '',
        kingdom = ifelse(is.na(kingdom) | kingdom == '', 'Unassigned', kingdom),
        specificEpithet = ifelse(!(is.na(Species_exact) | Species_exact == ''), Species_exact, specificEpithet),
        specificEpithet = ifelse( (!(is.na(genus) | genus == '')), str_replace(specificEpithet, paste('^',genus, '[_[:space:]]' ,sep=''), ''), specificEpithet),
        specificEpithet = ifelse( str_detect(specificEpithet, '^[sS]p{1,2}.?$'), '', specificEpithet),
        infraspecificEpithet = ifelse( str_detect(specificEpithet, '[_[:space:]]'), specificEpithet, ''),
        infraspecificEpithet = str_replace(infraspecificEpithet, '^[^_[:space:]]*[_[:space:]]', ''),
        specificEpithet = str_replace(specificEpithet, paste('[_[:space:]]', infraspecificEpithet ,sep=''), ''),
    ) %>%
    relocate(otu, .after = infraspecificEpithet) %>%
    relocate(associatedSequences, .before = domain) %>%
    select_if(!names(.) %in% c('confidence','domain', 'Species_exact', 'SH', 'BOLD_bin', 'Supergroup', 'Division', 'Subdivision')) %>%
    write_tsv("asv-table.tsv", na = '')
