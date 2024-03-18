#!/bin/sh

derep="$1"

# Untar any tar file in the working directory
tar xzf database.tar.gz

# Greengenes 13_8
if [ -d "gg_13_8_otus" ]; then
    mv gg_13_8_otus/rep_set/${derep}_otus.fasta gg_13_8_otus_rep_set_${derep}_otus.seq.fasta
    mv gg_13_8_otus/rep_set_aligned/${derep}_otus.fasta gg_13_8_otus_rep_set_aligned_${derep}_otus.alnseq.fasta
    mv gg_13_8_otus/taxonomy/${derep}_otu_taxonomy.txt gg_13_8_otus_taxonomy_${derep}_otu_taxonomy.tax.txt
    # remove uncompressed folder
    rm -r gg_13_8_otus
elif [ -d "SILVA_128_QIIME_release" ]; then
    mv SILVA_128_QIIME_release/rep_set/rep_set_all/${derep}/${derep}_otus.fasta SILVA_128_QIIME_release_rep_set_all_${derep}_otus.seq.fasta
    gunzip -c SILVA_128_QIIME_release/rep_set_aligned/${derep}/${derep}_otus_aligned.fasta.gz > SILVA_128_QIIME_release_rep_set_aligned_${derep}_otus_aligned.alnseq.fasta
    mv SILVA_128_QIIME_release/taxonomy/taxonomy_all/${derep}/consensus_taxonomy_7_levels.txt SILVA_128_QIIME_release_taxonomy_all_${derep}_consensus_taxonomy_7_levels.tax.txt
    # remove uncompressed folder
    rm -r SILVA_128_QIIME_release
else
    echo "No expected directory detected"
fi
