process BCFTOOLS_COMPRESS_INDEX {

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${vcf}.gz"), path("${vcf}.gz.tbi"), emit: indexed
    path "versions.yml",                                       emit: versions

    script:
    """
    bcftools view -O z $vcf > ${vcf}.gz

    bcftools index -t ${vcf}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${vcf}.gz
    touch ${vcf}.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

}