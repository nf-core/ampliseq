gzip -d *_toGenus.fa.gz 
cp *_toGenus.fa assignTaxonomy.fna

gzip -d *_toSpecies.fa.gz
cp *_toSpecies.fa addSpecies.fna
