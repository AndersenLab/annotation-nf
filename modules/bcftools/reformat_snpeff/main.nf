process BCFTOOLS_REFORMAT_SNPEFF {

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.species}.snpeff.tsv"), emit: tsv
    path "versions.yml",                                 emit: versions

    script:
    """
    # Create header for output file
    echo -e "CHROM\\tPOS\\tREF\\tALT\\tCONSEQUENCE\\tIMPACT\\tAA\\tSTRAIN\\tTRANSCRIPT" > ${meta.species}.snpeff.tsv

    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/ANN[\\t%SAMPLE=%GT]\\n' ${vcf} | \\
        awk -F'\\t' '{
            ALT_samples = "";  # Initialize string to collect samples with alt allele
            for (i = 6; i <= NF; i++) {  # Loop through fields containing sample=genotype
                if (\$i ~ /0\\/1|1\\/0|1\\/1/) {  # Check if sample has alt allele
                    sub(/=.*/, "", \$i);  # Remove genotype - leaving only sample name
                    if (ALT_samples == "") ALT_samples = \$i;
                    else ALT_samples = ALT_samples " " \$i;
                }
            }

            if (ALT_samples != "") {
                split(\$5, annotations, ",");  # Split multiple annotations into array
                for (j in annotations) {
                    split(annotations[j], snpEff, "|");  # Split each annotation by pipe

                    SnpEff_consequence = snpEff[2];
                    SnpEff_impact = snpEff[3];
                    SnpEff_AA_change = (snpEff[11] != "" ? snpEff[11] : "N/A");
                    transcript = snpEff[7]; #(snpEff[7] != "" ? snpEff[7] : "N/A");
                    if (SnpEff_consequence == "intergenic_region") {
                        transcript = "N/A";
                    }

                    print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" SnpEff_consequence "\\t" SnpEff_impact "\\t" SnpEff_AA_change "\\t" ALT_samples "\\t" transcript;
                }
            }
        }' >> ${meta.species}.snpeff.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.species}.snpeff.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

}