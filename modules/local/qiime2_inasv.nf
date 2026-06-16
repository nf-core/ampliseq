process QIIME2_INASV {
    tag "${asv}"
    label 'process_low'

    conda "${moduleDir}/envs/rachis-qiime2-linux-64-conda.yml"
    container "qiime2/qiime2:2026.4"

    input:
    path(asv)

    output:
    path("table.qza")    , emit: qza
    path "versions.yml"  , emit: versions_qiime2_inasv, topic: versions

    script:
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
