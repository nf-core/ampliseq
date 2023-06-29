process BARRNAP {
    tag "${fasta}"
    label 'process_low'

    conda "bioconda::barrnap=0.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/barrnap:0.9--hdfd78af_4' :
        'biocontainers/barrnap:0.9--hdfd78af_4' }"

    input:
    path(fasta)

    output:
    path( "*.matches.txt" ) , emit: matches
    path( "rrna.*.gff" )    , emit: gff
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def kingdom = task.ext.kingdom ?: "bac,arc,mito,euk"
    """
    IFS=',' read -r -a kingdom <<< \"$kingdom\"

    for KINGDOM in "\${kingdom[@]}"
    do
        barrnap \\
            --threads $task.cpus \\
            $args \\
            --kingdom \$KINGDOM \\
            --outseq Filtered.\${KINGDOM}.fasta \\
            < $fasta \\
            > rrna.\${KINGDOM}.gff

        #this fails when processing an empty file, so it requires a workaround!
        if [ -s Filtered.\${KINGDOM}.fasta ]; then
            grep -h '>' Filtered.\${KINGDOM}.fasta | sed 's/^>//' | sed 's/:\\+/\\t/g' | awk '{print \$2}' | sort -u >\${KINGDOM}.matches.txt
        else
            touch \${KINGDOM}.matches.txt
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        barrnap: \$(echo \$(barrnap --version 2>&1) | sed "s/^.*barrnap //g")
    END_VERSIONS
    """
}
