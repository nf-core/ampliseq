process SAVONT_ASV {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'

    conda "bioconda::savont=0.6.1c"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/savont:0.6.1c--hec9b1f2_0' :
        'biocontainers/savont:0.6.1c--hec9b1f2_0' }"

    input:
    tuple val(meta), path(reads), val(sample_string)

    output:
    path("*_feature-table.tsv") , emit: asv
    path("*_final_asvs.fasta")  , emit: fasta
    path("*_final_clusters.tsv"), emit: clusters
    path("*_stats.tsv")         , emit: stats
    path("savont_asv_*/temp/*") , emit: temp
    tuple val(meta), path("savont_asv_*"), emit: output_folder
    path("*_savont.log")        , emit: log
    tuple val("${task.process}"), val('savont'), eval("savont --version 2>&1 | cut -d ' ' -f 2"), topic: versions, emit: versions_savont

    script:
    def prefix = task.ext.prefix ?: "$meta.id"
    def args = task.ext.args ?: ''
    def reads_cmd = reads instanceof List ? "--pooled-samples ${reads.join(' ')}" : "${reads}"
    """
    savont asv \\
        -t $task.cpus \\
        $args \\
        -o savont_asv_${prefix} \\
        $reads_cmd

    # adjust the naming of the samples:
    sed '1s/.*/${sample_string}/' savont_asv_${prefix}/feature-table.tsv > ${prefix}_feature-table.tsv

    # output count per sample
    savont_asv_stats.sh ${prefix}_feature-table.tsv >${prefix}_stats.tsv

    # copy other files to include the prefix
    cp savont_asv_${prefix}/final_asvs.fasta ${prefix}_final_asvs.fasta
    cp savont_asv_${prefix}/final_clusters.tsv ${prefix}_final_clusters.tsv
    cp savont_asv_${prefix}/savont_*.log ${prefix}_savont.log
    """
}
