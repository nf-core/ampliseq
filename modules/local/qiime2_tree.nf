process QIIME2_TREE {
    label 'process_medium'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.8"

    input:
    path(repseq)

    output:
    path("rooted-tree.qza"), emit: qza
    path("tree.nwk")       , emit: nwk
    path "versions.yml"    , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    qiime alignment mafft \
        --i-sequences ${repseq} \
        --o-alignment aligned-rep-seqs.qza \
        --p-n-threads ${task.cpus}
    qiime alignment mask \
        --i-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza
    qiime phylogeny fasttree \
        --i-alignment masked-aligned-rep-seqs.qza \
        --p-n-threads ${task.cpus} \
        --o-tree unrooted-tree.qza
    qiime phylogeny midpoint-root \
        --i-tree unrooted-tree.qza \
        --o-rooted-tree rooted-tree.qza
    qiime tools export --input-path rooted-tree.qza  \
        --output-path phylogenetic_tree
    cp phylogenetic_tree/tree.nwk .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
    END_VERSIONS
    """
}
