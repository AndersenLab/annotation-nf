process LOCAL_GFF_SORT {

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${meta.species}.sorted_gff.tsv")

    script:
    """
    awk -F'\\t' '\$3 != "gene" && \$3 != "intron" && \$3 != "exon" && \$3 != "five_prime_UTR" && \$3 != "three_prime_UTR" && \$3 != "CDS" && \$3 != "start_codon" && \$3 != "stop_codon" {print \$9}' $gff | \
    awk -F'[;=]' '{
        transcript=""; wbgene="";gene=""
        for(i=1; i<=NF; i++) {
            if ($$i == "ID" && $$(i+1) ~ /^transcript:/) { transcript=$$(i+1); gsub("transcript:", "", transcript); }
            if ($$i == "ID" && $$(i+1) ~ /^Pseudogene:/) { transcript=$$(i+1); gsub("Pseudogene:", "", transcript); }
            if ($$i == "Parent" && $$(i+1) ~ /^gene:/) { wbgene=$$(i+1); gsub("gene:", "", wbgene); }
            if ($$i == "locus") { gene=$$(i+1); }
        }
        if (gene == "") { gene = "N/A"; }

        if (transcript != "" && wbgene != "") print transcript "\\t" wbgene "\\t" gene;
    }' | sort -t$'\t' -k1,1 > ${meta.species}.sorted_gff.tsv
    """

    stub:
    """
    touch ${meta.species}.sorted_gff.tsv
    """

}