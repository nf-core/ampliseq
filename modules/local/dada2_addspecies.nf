process DADA2_ADDSPECIES {
    tag "${taxtable},${database}"
    label 'process_high'
    label 'single_cpu'

    conda (params.enable_conda ? "bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    path(taxtable)
    path(database)
    val(outfile)

    output:
    path(outfile)       , emit: tsv
    path "versions.yml" , emit: versions
    path "*.args.txt"   , emit: args

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    set.seed(100) # Initialize random number generator for reproducibility

    taxtable <- readRDS(\"$taxtable\")

    tx <- addSpecies(taxtable, \"$database\", $args, verbose=TRUE)

    # Create a table with specified column order
    taxa <- data.frame(tx[,1:8])
    taxa\$Species <- tx[,'Species']
    taxa\$confidence <- tx[,'confidence']
    taxa\$sequence   <- taxtable\$sequence

    write.table(taxa, file = \"$outfile\", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    write.table('addSpecies\t$args', file = "addSpecies.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\\"${task.process}\\":", paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}
