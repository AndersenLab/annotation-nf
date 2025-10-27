process BCFTOOLS_CSQ_ANNOTATION {

    label "md_long"
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path genome
    tuple val(meta2), path(gff)

    output:
    tuple val(meta), path("${meta.species}.csq.vcf.gz"), path("${meta.species}.csq.vcf.gz.tbi"), emit: csq
    path "versions.yml",                                                                         emit: versions

    script:
    """
    bcftools csq \
        -O z --fasta-ref $genome \
        --gff-annot $gff \
        --ncsq 1000 \
        --phase a $vcf > ${meta.species}.csq.vcf.gz

    bcftools index -t ${meta.species}.csq.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.species}.csq.vcf.gz
    touch ${meta.species}.csq.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

}