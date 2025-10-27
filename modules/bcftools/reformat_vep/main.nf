process BCFTOOLS_REFORMAT_VEP {

    label "mid"

    input:
    tuple val(meta), path(vcf)
    val hdr_option

    output:
    tuple val(meta), path("${meta.species}.hdr${hdr_option}.vep.tsv"), emit: tsv
    path "versions.yml",                                               emit: versions

    script:
    """
    if [[ ${hdr_option} -eq 1 ]]; then
        HDR="YES"
    else
        HDR="NO"
    fi

    # Create header for output file
    echo -e "CHROM\\tPOS\\tREF\\tALT\\tCONSEQUENCE\\tIMPACT\\tAA\\tSTRAIN\\tHDR\\tTRANSCRIPT" > ${meta.species}.hdr${hdr_option}.vep.tsv

    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/CSQ[\\t%SAMPLE=%GT:%HDR]\\n' -i 'HDR=${hdr_option}' ${vcf} | \\
        awk -F'\\t' -v HDR=\$HDR '{
            ALT_samples = "";
            for (i = 6; i <= NF; i++) {
                if (\$i ~ /0\\/1|1\\/0|1\\/1/) {
                    sub(/=.*/, "", \$i);
                    if (ALT_samples == "") ALT_samples = \$i;
                    else ALT_samples = ALT_samples " " \$i;
                }
            }
            if (ALT_samples != "") {
                if (\$5 != ".") {
                    split(\$5, annotations, ",");
                    for (j in annotations) {
                        split(annotations[j], vep, "|");  #split VEP annotation by pipe
                        VEP_consequence = vep[2];
                        VEP_impact = vep[3];
                        VEP_AA_change = (vep[16] != "" ? vep[16] : "N/A");
                        transcript = (vep[7] != "" ? vep[7] : (vep[5] != "" ? vep[5] : "N/A"));
                        if (transcript ~ /&/) {
                            sub(/^.*&/, "", transcript);
                        }
                        # blosum62score = (vep[25] != "" ? vep[25] : "N/A");
                        print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" VEP_consequence "\\t" VEP_impact "\\t" VEP_AA_change "\\t" ALT_samples "\\t" HDR "\\t" transcript;
                    }
                } else {
                    print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\tN/A\\tN/A\\tN/A\\t" ALT_samples "\\t" HDR "\\tN/A";
                }
            }
        }' >> ${meta.species}.hdr${hdr_option}.vep.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.species}.hdr${hdr_option}.vep.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

}