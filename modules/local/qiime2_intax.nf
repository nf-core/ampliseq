process QIIME2_INTAX {
    tag "${tax}"
    label 'process_low'

    conda "${moduleDir}/envs/rachis-qiime2-linux-64-conda.yml"
    container "qiime2/qiime2:2026.4"

    input:
    path(tax) //ASV_tax_species.tsv
    val(script)

    output:
    path("taxonomy.qza") , emit: qza
    path "versions.yml"  , emit: versions_qiime2_intax, topic: versions

    script:
    def script_cmd = script ? "$script $tax" : "cp $tax tax.tsv"
    """
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    $script_cmd

    qiime tools import \\
        --type 'FeatureData[Taxonomy]' \\
        --input-format HeaderlessTSVTaxonomyFormat \\
        --input-path tax.tsv \\
        --output-path taxonomy.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
