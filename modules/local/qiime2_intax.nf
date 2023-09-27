process QIIME2_INTAX {
    tag "${tax}"
    label 'process_low'

    container "qiime2/core:2022.11"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path(tax) //ASV_tax_species.tsv
    val(script)

    output:
    path("taxonomy.qza") , emit: qza
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def script_cmd = script ? "$script $tax" : "cp $tax tax.tsv"
    """
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
