process NOVASEQ_ERR {
    tag "$meta.run"
    label 'process_medium'

    conda "bioconda::bioconductor-dada2=1.38.0 conda-forge::r-base=4.5.2 conda-forge::r-digest=0.6.39 conda-forge::tbb=2022.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81153df5d53322e6d91b2c4c9bc4da50774fb1d101ead002fe75bb75fc6f036c/data' :
        'community.wave.seqera.io/library/bioconductor-dada2_r-base_r-digest_tbb:38acac09bac46f36' }"

    input:
    tuple val(meta), path(errormodel)

    output:
    tuple val(meta), path("*.md.err.rds"), emit: errormodel
    tuple val(meta), path("*.md.err.pdf"), emit: pdf
    tuple val(meta), path("*.md.err.svg"), emit: svg
    tuple val(meta), path("*.md.err.convergence.txt"), emit: convergence
    path "versions.yml"                  , emit: versions

    script:
    if (!meta.single_end) {
        """
        novaseq_err_pe.r ${errormodel[0]} ${errormodel[1]} ${meta.run}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/\\s.*\$//')
        END_VERSIONS
        """
    } else {
        """
        novaseq_err_se.r ${errormodel} ${meta.run}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/\\s.*\$//')
        END_VERSIONS
        """
    }
}
