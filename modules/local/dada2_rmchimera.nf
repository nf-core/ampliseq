process DADA2_RMCHIMERA {
    tag "$meta.run"
    label 'process_medium'

    conda (params.enable_conda ? "bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    tuple val(meta), path(seqtab)

    output:
    tuple val(meta), path("*.ASVtable.rds"), emit: rds
    path "versions.yml"                    , emit: versions
    path "*.args.txt"                      , emit: args

    script:
    def args = task.ext.args ?: ''
    def no_samples    = meta.id.size()
    def first_sample  = meta.id.first()
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    seqtab = readRDS("${seqtab}")

    #remove chimera
    seqtab.nochim <- removeBimeraDenovo(seqtab, $args, multithread=$task.cpus, verbose=TRUE)
    if ( ${no_samples} == 1 ) { rownames(seqtab.nochim) <- "${first_sample}" }
    saveRDS(seqtab.nochim,"${meta.run}.ASVtable.rds")

    write.table('removeBimeraDenovo\t$args', file = "removeBimeraDenovo.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}
