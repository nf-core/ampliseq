process QIIME2_INASV {
    tag "${asv}"
    label 'process_low'

    container "qiime2/core:2023.7"

    input:
    path(asv)

    output:
    path("table.qza")    , emit: qza
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    # remove first line if needed
    sed '/^# Constructed from biom file/d' "$asv" > biom-table.txt

    # load into QIIME2
    biom convert -i biom-table.txt -o table.biom --table-type="OTU table" --to-hdf5
    qiime tools import \\
        --input-path table.biom \\
        --type 'FeatureTable[Frequency]' \\
        --input-format BIOMV210Format \\
        --output-path table.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
