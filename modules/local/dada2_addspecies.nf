// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process DADA2_ADDSPECIES {
    tag "${taxtable},${database}"
    label 'process_high'
    label 'single_cpu'
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
    path(taxtable)
    path(database)
    val(outfile)
    
    output:
    path(outfile)       , emit: tsv
    path "*.version.txt", emit: version
    path "*.args.txt"   , emit: args

    script:
    def software      = getSoftwareName(task.process)
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    set.seed(100) # Initialize random number generator for reproducibility

    taxtable <- readRDS(\"$taxtable\")

    tx <- addSpecies(taxtable, \"$database\", $options.args, verbose=TRUE)

    # Create a table with specified column order
    tmp <- data.frame(row.names(tx)) # To separate ASV_ID from sequence
    taxa <- data.frame(
        ASV_ID = tx[,"ASV_ID"],
        Kingdom = tx[,"Kingdom"],
        Phylum = tx[,"Phylum"],
        Class = tx[,"Class"],
        Order = tx[,"Order"],
        Family = tx[,"Family"],
        Genus = tx[,"Genus"],
        Species = tx[,"Species"],
        confidence = tx[,"confidence"],
	sequence = tmp[,],
	row.names=row.names(tmp)
    )

    write.table(taxa, file = \"$outfile\", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    write.table('addSpecies\t$options.args', file = "addSpecies.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    """
}
