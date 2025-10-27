process BCFTOOLS_HDR_MATRIX {

    label "md_mid"

    input:
    tuple val(meta), path(strain_matrix)
    path genos, stageAs: "genos/*"
    tuple val(meta1), path(vcf), path(vcf_tbi)

    output:
    tuple val(meta), path("HDR.vcf.gz"), path("HDR.vcf.gz.tbi"), emit: hdr
    path "versions.yml",                                         emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    # Get sample list from VCF
    bcftools query -l ${vcf} > sample_list.txt
    
    # Process each sample and add genotype data
    awk 'BEGIN{ printf "-d \\"\\\\t\\" ${strain_matrix}" }{ printf " genos/%s.geno", \$1 }' sample_list.txt | \\
        xargs paste > HDR_strain_matrix.tsv

    # Compress and index the matrix
    bgzip -c HDR_strain_matrix.tsv > HDR_strain_matrix.tsv.gz
    tabix -s1 -b2 -e2 HDR_strain_matrix.tsv.gz
    rm HDR_strain_matrix.tsv

    # Create HDR header for VCF
    echo '##FORMAT=<ID=HDR,Number=1,Type=Integer,Description="HDR annotation: 1 for YES, 0 for NO">' > HDR.header 

    # Add HDR annotations to the VCF
    bcftools annotate \\
        --annotations HDR_strain_matrix.tsv.gz \\
        --columns CHROM,POS,FORMAT/HDR \\
        --header-lines HDR.header \\
        --samples-file sample_list.txt \\
        --output HDR.vcf.gz \\
        --output-type z \\
        ${vcf}

    # Index the output VCF
    bcftools index -t HDR.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/^.*bcftools //')
        tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/^.*tabix //')
        bgzip: \$(bgzip --version 2>&1 | head -n1 | sed 's/^.*bgzip //')
    END_VERSIONS
    """

    stub:
    """
    touch HDR.vcf.gz
    touch HDR.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/^.*bcftools //')
        tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/^.*tabix //')
        bgzip: \$(bgzip --version 2>&1 | head -n1 | sed 's/^.*bgzip //')
    END_VERSIONS
    """
}