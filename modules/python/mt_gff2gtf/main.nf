process PYTHON_MT_GFF2GTF {

    input:
    tuple val(meta), path(gff)
    path "mt_gff2gtf.py"

    output:
    tuple val(meta), path("*.gtf"), emit: gtf
    path "versions.yml",            emit: versions

    script:
    def gtf = gff.name.lastIndexOf('.').with {it != -1 ? gff.name[0..<it] : gff.name} + "_mt.gtf"
    """
    python mt_gff2gtf.py ${gff} ${gtf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version |& sed '1!d; s/^.*Python //' )
    END_VERSIONS
    """

    stub:
    def gtf = gff.name.lastIndexOf('.').with {it != -1 ? gff.name[0..<it] : gff.name} + "_mt.gtf"
    """
    touch ${gtf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version |& sed '1!d; s/^.*Python //' )
    END_VERSIONS
    """
}