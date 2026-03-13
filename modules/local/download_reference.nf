process DOWNLOAD_REFERENCE {
    tag "${local_filename}"
    label 'process_single'

    conda "conda-forge::wget=1.21.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b54fa9135194c72a18d00db6b399c03248103f87e43ca75e4b50d61179994b3/data':
        'community.wave.seqera.io/library/wget:1.21.4--8b0fcde81c17be5e' }"

    input:
    val ref_url

    output:
    path "${local_filename}", emit: db
    // version output file is not wanted because of overwriting issues when storeDir is used

    script:
    local_filename = file(ref_url).name
    """
    wget -O "${local_filename}" "${ref_url}"
    """
}
