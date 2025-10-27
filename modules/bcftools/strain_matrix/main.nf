process BCFTOOLS_STRAIN_MATRIX {

    label "sm"

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("strain_matrix.tsv"),                      emit: matrix
    tuple val(meta), path("no_mt_vcf.bed"),                          emit: no_mt_bed
    tuple val(meta), path("no_mt.vcf.gz"), path("no_mt.vcf.gz.tbi"), emit: no_mt
    tuple val(meta), path("mt.vcf.gz"), path("mt.vcf.gz.tbi"),       emit: mt
    path "sample_names.txt",                                         emit: samples
    path "versions.yml",                                             emit: versions
    
    shell:
    """
    # Filtering VCF to only biallelic sites
    bcftools view -m2 -M2 -e 'CHROM=="MtDNA"' -O z -o no_mt.vcf.gz $vcf
    bcftools index -t no_mt.vcf.gz
    
    bcftools view -m2 -M2 -i 'CHROM=="MtDNA"' -O z -o mt.vcf.gz $vcf 
    bcftools index -t mt.vcf.gz

    # Extract list of strain names
    sample_names=(\$(bcftools query -l $vcf))
    for name in \${sample_names[*]}; do
        echo \${name} >> sample_names.txt
    done

    # Convert VCF to BED format
    bcftools query -f '%CHROM\\t%POS\\n' no_mt.vcf.gz | awk '{print \$1"\\t"\$2"\\t"\$2+1}' > no_mt_vcf.bed

    # Create output file genotype matrix
    cut -f1,2 no_mt_vcf.bed > strain_matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch no_mt.vcf.gz
    touch no_mt.vcf.gz.tbi
    touch mt.vcf.gz
    touch mt.vcf.gz.tbi
    touch strain_matrix.tsv
    touch sample_names.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}