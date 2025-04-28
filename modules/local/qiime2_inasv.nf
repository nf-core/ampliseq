process QIIME2_INASV {
    tag "${asv}"
    label 'process_low'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    path(asv)

    output:
    path("table.qza")    , emit: qza
    path "versions.yml"  , emit: versions

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
