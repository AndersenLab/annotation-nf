process BCFTOOLS_REFORMAT_ANNOVAR {

    label "mid"

    input:
    tuple val(meta), path(vcf)
    val hdr_option

    output:
    tuple val(meta), path("${meta.species}.hdr${hdr_option}.annovar.tsv"), emit: tsv
    path "versions.yml",                                                   emit: versions

    script:
    """
    if [[ ${hdr_option} -eq 1 ]]; then
        HDR="YES"
    else
        HDR="NO"
    fi

    # Create header for output file
    echo -e "CHROM\\tPOS\\tREF\\tALT\\tCONSEQUENCE\\tIMPACT\\tAA\\tSTRAIN\\tHDR\\tTRANSCRIPT" > ${meta.species}.hdr${hdr_option}.annovar.tsv

    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO[\\t%SAMPLE=%GT:%HDR]\\n' -i 'HDR=${hdr_option}' ${vcf} | \\
        awk -F'\\t' -v HDR=\$HDR '{
            # Extract samples with alternative alleles
            ALT_samples = "";
            for (i = 6; i <= NF; i++) {
                if (\$i ~ /0\\/1|1\\/0|1\\/1/) {
                    sample_name = \$i;
                    sub(/=.*/, "", sample_name);
                    if (ALT_samples == "") ALT_samples = sample_name;
                    else ALT_samples = ALT_samples " " sample_name;
                }
            }
            
            # Only process variants with alternative samples
            if (ALT_samples != "") {
                # Initialize annotation variables
                ANV_consequence = "N/A";
                ANV_impact = "N/A";
                ANV_AA_change = "N/A";
                transcript = "N/A";
                transcripts = "";
                genes = "";
                reported = 0;
                
                # Parse INFO field annotations
                if (\$5 != "" && \$5 != ".") {
                    split(\$5, anv, ";");
                    for (i = 1; i <= length(anv); i++) {
                        if (index(anv[i], "=") > 0) {
                            split(anv[i], fields, "=");
                            tag = fields[1];
                            value = (length(fields) > 1) ? fields[2] : "";
                            
                            if (tag == "Func.refGene" && value != "") {
                                ANV_consequence = (value == "." ? "N/A" : value);
                            }
                            if (tag == "ExonicFunc.refGene" && value != "") {
                                ANV_impact = (value == "." ? "N/A" : value);
                            }
                            if (tag == "AAChange.refGene" && value != "" && value != ".") {
                                transcripts = value;
                            }
                            if (tag == "Gene.refGene" && value != "" && value != ".") {
                                genes = value;
                            }
                        }
                    }
                }
                if (length(transcripts) > 0) {
                    # Handle multiple transcripts separated by commas
                    split(transcripts, transcript_list, "\\\\\\\\x3b");
                    for (j = 1; j <= length(transcript_list); j++) {
                        reported = 1;
                        if (transcript_list[j] != "") {
                            split(transcript_list[j], AA_change_parts, ":");
                            current_AA_change = "N/A";
                            current_transcript = "N/A";
                            
                            for (k = 1; k <= length(AA_change_parts); k++) {
                                # Extract amino acid change (p.XXX)
                                if (AA_change_parts[k] ~ /^p\\./) current_AA_change = AA_change_parts[k];
                                if (AA_change_parts[k] == "transcript") current_transcript = AA_change_parts[k+1];
                            }
                            
                            # Print one line per transcript
                            print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" ANV_consequence "\\t" ANV_impact "\\t" current_AA_change "\\t" ALT_samples "\\t" HDR "\\t" current_transcript;
                        }
                    }
                } else if (length(genes) > 0) {
                    # Handle multiple transcripts separated by commas
                    split(genes, gene_list, "\\\\\\\\x3b");
                    split(ANV_consequence, consequence_list, "\\\\\\\\x3b");
                    for (j = 1; j <= length(gene_list); j++) {
                        reported = 1;
                        if (gene_list[j] != "") {
                            split(gene_list[j], gene_parts, ":");
                            current_gene = gene_parts[2];
                            if (j > length(consequence_list)) ANV_consequence = consequence_list[1];
                            else ANV_consequence = consequence_list[j];

                            # Print one line per gene
                            print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" ANV_consequence "\\t" ANV_impact "\\t" ANV_AA_change "\\t" ALT_samples "\\t" HDR "\\t" current_gene;
                        }
                    }
                }
                
                # If no AAChange.refGene field found, still output the variant
                if (reported == 0) {
                    print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"ANV_consequence"\\t"ANV_impact"\\t"ANV_AA_change"\\t"ALT_samples"\\t"HDR"\\t"transcript;
                }
            }
        }' >> ${meta.species}.hdr${hdr_option}.annovar.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    echo -e "CHROM\\tPOS\\tREF\\tALT\\tFUNC\\tEXONIC_FUNC\\tAA_CHANGE\\tALT_SAMPLES\\tHDR\\tTRANSCRIPT" > ${meta.species}.hdr${hdr_option}.annovar.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

}