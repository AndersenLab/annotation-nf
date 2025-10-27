process LOCAL_GFF2TSV {

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("gff.tsv")
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Extract relevant features from GFF and parse attributes
    cat ${gff} | \\
    awk -F'\\t' '
        \$3 != "gene" && 
        \$3 != "intron" && 
        \$3 != "exon" && 
        \$3 != "five_prime_UTR" && 
        \$3 != "three_prime_UTR" && 
        \$3 != "CDS" && 
        \$3 != "start_codon" && 
        \$3 != "stop_codon" &&
        \$3 != "tRNA" {
            print \$9
        }' | \\
    awk -F'[;=]' '
        {
            transcript = ""
            wbgene = ""
            gene = ""
            
            # Parse the attributes field
            for(i = 1; i <= NF; i++) {
                if (\$i == "ID" && \$(i+1) ~ /^transcript:/) {
                    transcript = \$(i+1)
                    gsub("transcript:", "", transcript)
                }
                if (\$i == "ID" && \$(i+1) ~ /^Pseudogene:/) {
                    transcript = \$(i+1)
                    gsub("Pseudogene:", "", transcript)
                }
                if (\$i == "Parent" && \$(i+1) ~ /^gene:/) {
                    wbgene = \$(i+1)
                    gsub("gene:", "", wbgene)
                }
                if (\$i == "locus") {
                    gene = \$(i+1)
                }
            }
            
            # Set default value for gene if empty
            if (gene == "") {
                gene = "N/A"
            }
            
            # Output if we have both transcript and wbgene
            if (transcript != "" && wbgene != "") {
                print transcript "\\t" wbgene "\\t" gene
            }
        }
        END {
            print "N/A\\tN/A\\tN/A"
        }' > unsorted.tsv
    
    # Sort the output by first column
    sort -t\$'\\t' -k1,1 unsorted.tsv > gff.tsv
    
    # Clean up temporary file
    rm -f unsorted.tsv
    """

    stub:
    """
    touch gff.tsv
    """
}