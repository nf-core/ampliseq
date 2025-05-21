#!/bin/bash

# The file for taxonomy assignment, identified by containing "toSpecies" in the name
gunzip -c *toSpecies*gz > assignTaxonomy.fna

# and the file for add species, identified by containing "assign" in the name, is renamed
mv *assign*gz addSpecies.fna.gz
