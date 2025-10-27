process JAVA_SNPEFF_ANNOTATION {

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path genome
    tuple val(meta2), path(gtf)
    path "base.config"

    output:
    tuple val(meta), path("${meta.species}.snpeff.vcf"), emit: snpeff
    path "versions.yml",                                 emit: versions

    script:
    """
    cp base.config snpEff.config

    GENOME_NAME=${meta.species}
    echo "\${GENOME_NAME}.genome : ${meta.species}" >> snpEff.config
    echo "\${GENOME_NAME}.MtDNA.codonTable = Invertebrate_Mitochondrial" >> snpEff.config

    mkdir ${meta.species}
    ln -s \${PWD}/${genome} ${meta.species}/sequences.fa
    ln -s \${PWD}/${gtf} ${meta.species}/genes.gtf

    snpEff build -noCheckCds -noCheckProtein -gtf22 -v \${GENOME_NAME}

    snpEff eff -csvStats snpeff.stats.csv \
              -nodownload \
              -dataDir . \
              -config snpEff.config \
              \${GENOME_NAME} \
              ${vcf} > ${meta.species}.snpeff.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$( java --version |& head -n 1 |& cut -f2 )
    "${task.process}":
        snpEff: \$( snpEff  -version |& cut -f2 )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.species}.snpeff.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$( java --version |& head -n 1 |& cut -f2 )
    "${task.process}":
        snpEff: \$( snpEff  -version |& cut -f2 )
    END_VERSIONS
    """

}