#!/bin/bash

for f in $(ls);
do
    c=$(gunzip -c $f | head -1 | wc -m | egrep -o "[0-9]+")
    echo -e "$f\t$c" >> tmp
done

assignTaxonomy=$(cat tmp | sort -k 2 -n -r | head -1 | cut -f1)
gunzip -c $assignTaxonomy > assignTaxonomy.fna

addSpecies=$(cat tmp | sort -k 2 -n | head -1 | cut -f1)
gunzip -c $addSpecies > addSpecies.fna
