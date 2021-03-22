// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process DADA2_TAXONOMY {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconductor-dada2=1.18.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.18.0--r40h5f743cb_0"
    } else {
        container "quay.io/biocontainers/bioconductor-dada2:1.18.0--r40h5f743cb_0"
    }

    input:
    path(fasta)
    path(database)
    
    output:
    path( "ASV_tax.tsv" ), emit: tsv
    path( "ASV_tax.rds" ), emit: rds
    path "*.version.txt" , emit: version
    path "*.args.txt"    , emit: args

    script:
    def software      = getSoftwareName(task.process)
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    set.seed(100) # Initialize random number generator for reproducibility

    seq <- getSequences(\"$fasta\", collapse = TRUE, silence = FALSE)
    taxa <- assignTaxonomy(seq, \"$database\", $options.args, multithread = $task.cpus, verbose=TRUE)
    saveRDS(taxa, "ASV_tax.rds")

    taxa <- cbind(sequence = rownames(taxa), taxa)
    write.table(taxa, file = "ASV_tax.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    write.table('assignTaxonomy\t$options.args', file = "assignTaxonomy.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    """
}