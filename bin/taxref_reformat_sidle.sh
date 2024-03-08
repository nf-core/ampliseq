#!/bin/sh

# Untar any tar file in the working directory
tar xzf database.tar.gz

# Greengenes 13_8
if [ -d "gg_13_8_otus" ]; then
    mv gg_13_8_otus/rep_set/99_otus.fasta gg_13_8_otus_rep_set_99_otus.seq.fasta
    mv gg_13_8_otus/rep_set_aligned/99_otus.fasta gg_13_8_otus_rep_set_aligned_99_otus.alnseq.fasta
    mv gg_13_8_otus/taxonomy/99_otu_taxonomy.txt gg_13_8_otus_taxonomy_99_otu_taxonomy.tax.txt
    # remove uncompressed folder
    rm -r gg_13_8_otus
elif [ -d "SILVA_128_QIIME_release" ]; then
    mv SILVA_128_QIIME_release/rep_set/rep_set_all/99/99_otus.fasta SILVA_128_QIIME_release_rep_set_all_99_otus.seq.fasta
    gunzip -c /SILVA_128_QIIME_release/rep_set_aligned/99/99_otus_aligned.fasta.gz > SILVA_128_QIIME_release_rep_set_aligned_99_otus_aligned.alnseq.fasta
    mv SILVA_128_QIIME_release/taxonomy/taxonomy_all/99/consensus_taxonomy_7_levels SILVA_128_QIIME_release_taxonomy_all_99_consensus_taxonomy_7_levels.tax.txt
    # remove uncompressed folder
    rm -r SILVA_128_QIIME_release
else
    echo "No expected directory detected"
fi





