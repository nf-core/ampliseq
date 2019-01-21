# nf-core/ampliseq: Troubleshooting

## Input files not found

If no file or less files than expected are picked up then something is wrong with your input file declaration (path) or with the naming of the files.

1. The path must be enclosed in quotes (`'` or `"`)
2. The files have to be or mimic Casava 1.8 paired-end demultiplexed fastq files with the naming sheme "[a-zA-Z0-9-]+_[a-zA-Z0-9-]+_L[0-9][0-9][0-9]_R{1,2}_001.fastq.gz". This is currently a limitation of QIIME2 file import.

If the pipeline can't find your files then you will get the following error

```
ERROR ~ Cannot find any reads matching: "[folder]/*_L[0-9][0-9][0-9]_R{1,2}_001.fastq.gz"
```


## Data organization
The pipeline can't take a list of multiple input files - it takes a single folder and picks up all files that match the pattern: 

`[a-zA-Z0-9-]+_[a-zA-Z0-9-]+_L[0-9][0-9][0-9]_R{1,2}_001.fastq.gz`

 If the input files do not follow the naming scheme, a directory with symlinks named as required linking to your actual data might be a solution. If your input files are scattered in different paths then we recommend that you generate a directory with symlinked files.

## Required computational resources

For learning a classifier from scratch for the appropriate primer pair (default) the pipeline requires 35GB memory. Other pipeline steps require a variable amount of memory increasing with number of samples, number of sequencing reads, number of sequencing runs, and sample complexity (number of unique sequences).

When memory is limiting you can use a pre-trained classifier with `--classifier [Path/To/Classifier.qza]` and e.g. `--max_memory [15.GB]`.

When number of CPUs are a matter of concern use e.g.`--max_cpus 4`.

Detailed decribtions can be found in [Running the pipeline](usage.md)

## Error related to metadata sheet

Please have a look at [QIIME2 metadata requirements](https://docs.qiime2.org/2018.6/tutorials/metadata). 
Generally, do not use whitespace characters or any special character.
The file has to contain tab-separated values. The first column should be named "ID" and contain unique terms which represent the beginning (usually referred to as sample id) of the sequencing files. 
For example the files L2S357_15_L001_R1_001.fastq.gz and L2S357_15_L001_R2_001.fastq.gz would have the ID L2S357.
Column names have to be unique and cannot be empty. Metadata values can be categorical (text) or numeric (containing only numbers or are empty). Empty cells represent missing data.

## Low read count after DADA2
This can have several reasons:
1. Poor data quality: High lost at quality filtering
2. Merging inefficient: Too less overlap for merging
3. Unusual high percentage of chimeric sequences: Remaining primer sequences

Solutions might be (numbered as above):
1. More agressive truncation values
2. Allow more overlap using higher truncation values
3. Double check the primer sequence

## "ValueError: CategoricalMetadataColumn does not support strings with leading or trailing whitespace characters"
This is a shortcoming of the QIIME2 classifier trainer and QIIME2 to handle '#' when specific sequences are in a sample related to specific taxonomies of the classifier. 

The solution is to remove all hash signs from the taxonomy strings using:

`--classifier_removeHash`

Also see [Running the pipeline](usage.md).

## Extra resources and getting help
If you still have an issue with running the pipeline then feel free to contact us.
Have a look at the [pipeline website](https://github.com/nf-core/ampliseq) to find out how.

If you have trouble interpreting results or with the methods itself, the [QIIME2 forum](https://forum.qiime2.org/) provides excellent help. Be aware that QIIME2 is in active development and you might use an older version of QIIME2 right now. 

If you have problems that are related to Nextflow and not our pipeline then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).
