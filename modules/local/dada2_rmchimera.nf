process DADA2_RMCHIMERA {
    tag "$meta.run"
    label 'process_medium'

    conda "bioconda::bioconductor-dada2=1.38.0 conda-forge::r-base=4.5.2 conda-forge::r-digest=0.6.39 conda-forge::tbb=2022.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81153df5d53322e6d91b2c4c9bc4da50774fb1d101ead002fe75bb75fc6f036c/data' :
        'community.wave.seqera.io/library/bioconductor-dada2_r-base_r-digest_tbb:38acac09bac46f36' }"

    input:
    tuple val(meta), path(seqtab)

    output:
    tuple val(meta), path("*.ASVtable.rds"), emit: rds
    path "versions.yml"                    , emit: versions
    path "*.args.txt"                      , emit: args

    script:
    def prefix = task.ext.prefix ?: "prefix"
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
    saveRDS(seqtab.nochim,"${prefix}.ASVtable.rds")

    write.table('removeBimeraDenovo\t$args', file = "removeBimeraDenovo.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}
