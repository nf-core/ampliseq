process QIIME2_INTAX {
    tag "${tax}"
    label 'process_low'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    path(tax) //ASV_tax_species.tsv
    val(script)

    output:
    path("taxonomy.qza") , emit: qza
    path "versions.yml"  , emit: versions

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
