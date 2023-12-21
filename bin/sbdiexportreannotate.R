#!/usr/bin/env Rscript

# sbdiexportreannotate.R
#
# A script that collates taxonomy data from Ampliseq to produce an updated
# annotation tsv file as close to ready for submission to the Swedish
# Biodiversity Data Infrastructure (SBDI) as possible.
#
# The script expects the following arguments: dbversion, ASV_tax_species.tsv, wfversion, predfile
#
# Author: daniel.lundin@lnu.se

suppressPackageStartupMessages(library(tidyverse))

# Get dbversion and taxonomy file from the command line
args            <- commandArgs(trailingOnly=TRUE)
dbversion       <- args[1]
taxfile         <- args[2]
taxmethod       <- args[3]
wfversion       <- args[4]
cut_its         <- args[5]
predfile        <- args[6]

cut_its = ifelse(cut_its == 'none', '', paste(' cut_its:', cut_its, sep=''))

# Read taxonomy table
taxonomy <- read.delim(taxfile, sep = '\t', stringsAsFactors = FALSE)

# Read the predictions table if provided, otherwise create one
if ( ! is.na(predfile) ) {
    predictions <- read.delim(predfile, sep = '\t', stringsAsFactors = FALSE)
    if ( nrow(predictions) < 1 ) {
        colnames <- names(predictions)
        predictions <- data.frame(ASV_ID = taxonomy$ASV_ID)
        predictions[,colnames[colnames != "ASV_ID"]] <- NA
    }
} else {
    predictions <- data.frame(ASV_ID = taxonomy$ASV_ID)
}

# Make sure it's congruent with the taxonomy table
predictions <- data.frame(
    'ASV_ID' = taxonomy$ASV_ID
) %>%
    left_join(predictions, by = 'ASV_ID' ) %>%
    distinct(ASV_ID, .keep_all = TRUE) %>%
    arrange(ASV_ID)

# Join tables and create missing columns
taxtable  <- taxonomy %>%
    inner_join(predictions, by = 'ASV_ID') %>%
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
    relocate(otu, .after = Species_exact) %>%
    rename_with(tolower, Domain:Species_exact) %>%
    rename(
        asv_id_alias = ASV_ID,
        asv_sequence = sequence,
        specificEpithet = species,
        annotation_confidence = confidence
    ) %>%
    mutate(across(.fns = ~str_replace_all(.,' ','_'))) %>%
    mutate(
        specificEpithet = ifelse(!(is.na(species_exact) | species_exact == ''), species_exact, specificEpithet),
        specificEpithet = ifelse( (!(is.na(genus) | genus == '')), str_replace(specificEpithet, paste('^',genus, '[_[:space:]]' ,sep=''), ''), specificEpithet),
        specificEpithet = ifelse( str_detect(specificEpithet, '^[sS]p{1,2}.?$'), '', specificEpithet),
        annotation_confidence =  ifelse((is.na(annotation_confidence) | annotation_confidence == ''), 0, annotation_confidence),
        scientificName = case_when(
            !(is.na(otu) | otu == '')                         ~ sprintf("%s", otu),
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
            !(is.na(otu) | otu == '')                         ~ 'unranked',
            !(is.na(specificEpithet) | specificEpithet == '') ~ 'species',
            !(is.na(genus)   | genus == '')                   ~ 'genus',
            !(is.na(family)  | family == '')                  ~ 'family',
            !(is.na(order)   | order == '')                   ~ 'order',
            !(is.na(class)   | class == '')                   ~ 'class',
            !(is.na(phylum)  | phylum == '')                  ~ 'phylum',
            !(is.na(kingdom) | kingdom == '')                 ~ 'kingdom',
            TRUE                                              ~ 'kingdom'
        ),
        domain = str_remove(domain, 'Reversed:_'),
        infraspecificEpithet = ifelse( str_detect(specificEpithet, '[_[:space:]]'), specificEpithet, ''),
        infraspecificEpithet = str_replace(infraspecificEpithet, '^[^_[:space:]]*[_[:space:]]', ''),
        specificEpithet = str_replace(specificEpithet, paste('[_[:space:]]', infraspecificEpithet ,sep=''), ''),
        date_identified = as.character(lubridate::today()),
        reference_db = dbversion,
        annotation_algorithm = case_when(
            (taxmethod == 'sintax')                         ~ paste('Ampliseq ',wfversion,' (https://nf-co.re/ampliseq) VSEARCH:sintax',cut_its, sep=' '),
            (!(is.na(otu) | otu == ''))                     ~ paste('Ampliseq ',wfversion,' (https://nf-co.re/ampliseq) addsh',cut_its, sep=' '),
            (!(is.na(species_exact) | species_exact == '')) ~ paste('Ampliseq ',wfversion,' (https://nf-co.re/ampliseq) DADA2:assignTaxonomy:addSpecies',cut_its, sep=' '),
            TRUE                                            ~ paste('Ampliseq ',wfversion,' (https://nf-co.re/ampliseq) DADA2:assignTaxonomy',cut_its, sep='')
        ),
        identification_references = 'https://docs.biodiversitydata.se/analyse-data/molecular-tools/#taxonomy-annotation',
        taxon_remarks = ifelse(!(is.na(domain) | domain == ''), paste('Domain = \'',domain,'\'',sep=''),''),
        kingdom = ifelse(is.na(kingdom) | kingdom == '', 'Unassigned', kingdom)
    ) %>%
    relocate(asv_sequence, .after = asv_id_alias) %>%
    relocate(scientificName:taxonRank, .after = asv_sequence) %>%
    relocate(infraspecificEpithet, .after = specificEpithet) %>%
    relocate(annotation_confidence, .after = otu) %>%
    relocate(date_identified:taxon_remarks, .after = annotation_confidence) %>%
    select_if(!names(.) %in% c('domain', 'species_exact', 'SH', 'BOLD_bin', 'Supergroup', 'Division', 'Subdivision')) %>%
    write_tsv("annotation.tsv", na = '')
