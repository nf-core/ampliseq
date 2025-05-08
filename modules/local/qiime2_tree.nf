process QIIME2_TREE {
    label 'process_medium'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    path(repseq)

    output:
    path("rooted-tree.qza"), emit: qza
    path("tree.nwk")       , emit: nwk
    path "versions.yml"    , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime alignment mafft \\
        --i-sequences ${repseq} \\
        --o-alignment aligned-rep-seqs.qza \\
        --p-n-threads ${task.cpus}
    qiime alignment mask \\
        --i-alignment aligned-rep-seqs.qza \\
        --o-masked-alignment masked-aligned-rep-seqs.qza
    qiime phylogeny fasttree \\
        --i-alignment masked-aligned-rep-seqs.qza \\
        --p-n-threads ${task.cpus} \\
        --o-tree unrooted-tree.qza
    qiime phylogeny midpoint-root \\
        --i-tree unrooted-tree.qza \\
        --o-rooted-tree rooted-tree.qza
    qiime tools export \\
        --input-path rooted-tree.qza  \\
        --output-path phylogenetic_tree
    cp phylogenetic_tree/tree.nwk .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
