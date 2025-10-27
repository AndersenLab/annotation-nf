process BCFTOOLS_REFORMAT_CSQ {

    label "mid"

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    val hdr_option

    output:
    tuple val(meta), path("${meta.species}.hdr${hdr_option}.csq.tsv"), emit: tsv
    path "versions.yml",                                               emit: versions

    script:
    """
    if [[ ${hdr_option} -eq 1 ]]; then
        HDR="YES"
    else
        HDR="NO"
    fi

    # Create header for output file
    echo -e "CHROM\\tPOS\\tREF\\tALT\\tCONSEQUENCE\\tAA\\tDNAchange\\tSTRAIN\\tHDR\\tTRANSCRIPT" > ${meta.species}.hdr${hdr_option}.csq.tsv

    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/BCSQ[\\t%SAMPLE=%GT:%HDR]\\n' -i 'HDR=${hdr_option}' ${vcf} | \\
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
                # Grantham_score = (\$5 != "." && \$5 != "" ? \$5 : "N/A");  # Default from column 5
                if (\$5 != ".") {
                    split(\$5, csq_entries, ",");  # Split multiple annotations by ","
                    for (j in csq_entries) {
                        split(csq_entries[j], csq, "|");  # Extract annotation components
                        # temp_Grantham = Grantham_score;  # Default to original score
                        if (csq[1] ~ /^@/) {
                            CSQ_consequence = csq[1];
                            CSQ_AA_change = csq[1];
                            DNA_change = csq[1];
                            # temp_Grantham = csq[1];  # Override only for this row
                            transcript = csq[1];
                        } else {
                            CSQ_consequence = (csq[1] != "" ? csq[1] : "N/A");
                            CSQ_AA_change = (csq[6] != "" ? csq[6] : "N/A");
                            DNA_change = (csq[7] != "" ? csq[7] : "N/A");
                            transcript = (csq[3] != "" ? csq[3] : (csq[2] != "" ? csq[2] : "N/A"));
                        }
                        print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"CSQ_consequence"\\t"CSQ_AA_change"\\t"DNA_change"\\t"ALT_samples"\\t"HDR"\\t"transcript;
                    }
                } else {
                    print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\tN/A\\tN/A\\tN/A\\t"ALT_samples"\\t"HDR"\\tN/A";
                }
            }
        }' >> ${meta.species}.hdr${hdr_option}.csq.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.species}.hdr${hdr_option}.csq.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

}