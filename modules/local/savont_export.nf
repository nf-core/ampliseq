process SAVONT_EXPORT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::savont=0.6.3"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/savont:0.6.3--hec9b1f2_0' :
        'biocontainers/savont:0.6.3--hec9b1f2_0' }"

    input:
    tuple val(meta), path(output_folders), val(sample_string)

    output:
    path("*_feature-table.tsv") , emit: asv
    path("*_final_asvs.fasta")  , emit: fasta
    path("*_savont_export.log") , emit: log
    tuple val("${task.process}"), val('savont'), eval("savont --version 2>&1 | cut -d ' ' -f 2"), topic: versions, emit: versions_savont

    script:
    def prefix = task.ext.prefix ?: "$meta.id"
    def args = task.ext.args ?: ''
    def folders = "${output_folders.join(' ')}"
    """
    savont export \\
        $args \\
        -i $folders \\
        -o savont_export \\
        --relabel $sample_string

    # copy other files to include the prefix
    cp savont_export/merged_feature_table.tsv ${prefix}_feature-table.tsv
    cp savont_export/merged_rep_seqs.fasta ${prefix}_final_asvs.fasta
    cp savont_export/savont_*.log ${prefix}_savont_export.log
    """
}
