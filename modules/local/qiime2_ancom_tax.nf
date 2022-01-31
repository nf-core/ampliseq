process QIIME2_ANCOM_TAX {
    tag "${table.baseName} - taxonomic level: ${taxlevel}"
    label 'process_medium'
    label 'single_cpu'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.8"

    input:
    tuple path(metadata), path(table), path(taxonomy) ,val(taxlevel)

    output:
    path "ancom/*"      , emit: ancom
    path "versions.yml" , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"
    mkdir ancom

    # Sum data at the specified level
    qiime taxa collapse \
            --i-table ${table} \
            --i-taxonomy ${taxonomy} \
            --p-level ${taxlevel} \
            --o-collapsed-table lvl${taxlevel}-${table}

    # Extract summarised table and output a file with the number of taxa
    qiime tools export --input-path lvl${taxlevel}-${table} --output-path exported/
    biom convert -i exported/feature-table.biom -o ${table.baseName}-level-${taxlevel}.feature-table.tsv --to-tsv

    if [ \$(grep -v '^#' -c ${table.baseName}-level-${taxlevel}.feature-table.tsv) -lt 2 ]; then
        echo ${taxlevel} > ancom/\"WARNING Summing your data at taxonomic level ${taxlevel} produced less than two rows (taxa), ANCOM can't proceed -- did you specify a bad reference taxonomy?\".txt
    else
        qiime composition add-pseudocount \
                --i-table lvl${taxlevel}-${table} \
                --o-composition-table comp-lvl${taxlevel}-${table}
        qiime composition ancom \
                --i-table comp-lvl${taxlevel}-${table} \
                --m-metadata-file ${metadata} \
                --m-metadata-column ${table.baseName} \
                --o-visualization comp-lvl${taxlevel}-${table.baseName}.qzv
        qiime tools export --input-path comp-lvl${taxlevel}-${table.baseName}.qzv \
                --output-path ancom/Category-${table.baseName}-level-${taxlevel}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
    END_VERSIONS
    """
}
