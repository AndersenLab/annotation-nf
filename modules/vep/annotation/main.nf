process VEP_ANNOTATION {

    label "md_long"

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path genome
    tuple val(meta2), path(gff)

    output:
    tuple val(meta), path("${meta.species}.vep.vcf"), emit: vep
    path "versions.yml",                              emit: versions

    script:
    """
    bgzip -c ${gff} > sample.gff.gz
    tabix -p gff sample.gff.gz

    vep \
    --vcf \
    --force_overwrite \
    --input_file $vcf \
    --fasta $genome \
    --gff sample.gff.gz \
    --output_file ${meta.species}.vep.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$( vep |& grep ensembl-vep |& awk '{print \$3}' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.species}.vep.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$( vep |& grep ensembl-vep |& awk '{print \$3}' )
    END_VERSIONS
    """

}