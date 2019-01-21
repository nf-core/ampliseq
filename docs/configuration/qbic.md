# nfcore/ampliseq: QBiC Configuration


QBiC has two clusters available to run and analyze data directly. For both systems, an appropriate configuration profile is existing, enabling a direct usability of the pipeline on the respective cluster environment. Until now, BINAC is the routinely used system.

# BINAC

You may use the pipeline with the `-profile binac` switch when starting the pipeline. A typical call could work like this for example
```
nextflow run nf-core/ampliseq \
    -profile binac,singularity \
    --reads "data" \
    --FW_primer GTGYCAGCMGCCGCGGTAA \
    --RV_primer GGACTACNVGGGTWTCTAAT \
    --metadata "data/Metadata.tsv"
``` 

---

[![QBiC](images/QBiC_logo.png)](https://portal.qbic.uni-tuebingen.de/portal/)

---
