process PERL_ANNOVAR_ANNOTATION {

    label "mem_long"

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path genome
    tuple val(meta2), path(gff)

    output:
    tuple val(meta), path("${meta.species}.annovar.vcf"), emit: annovar
    path "versions.yml",                                  emit: versions

    script:
    """
    mkdir db_build

    # convert gff to gtf
    gffread $gff -T -o db_build/${meta.species}.gtf

    # convert gtf format for Annovar
    gtfToGenePred -genePredExt db_build/${meta.species}.gtf db_build/${meta.species}_refGene.txt

    # create database from references
    retrieve_seq_from_fasta.pl --format refGene --seqfile $genome db_build/${meta.species}_refGene.txt --out db_build/${meta.species}_refGeneMrna.fa

    # create annovar annotations
    table_annovar.pl \
        $vcf db_build \
        --buildver ${meta.species} \
        --protocol refGene \
        --remove \
        --outfile annovar \
        --operation g \
        --thread ${task.cpus - 1} \
        --gff3dbfile $gff \
        --nopolish \
        --vcfinput

    mv annovar.${meta.species}_multianno.vcf ${meta.species}.annovar.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$( perl -e "print \$^V;" )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.species}.annovar.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$( perl -e "print \$^V;" )
    END_VERSIONS
    """

}