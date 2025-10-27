process LOCAL_CSQ_GFF {

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${meta.species}.csq.gff")

    script:
    """
    awk '{ \
        if ($3 == "mRNA"){ \
            gsub("product=", "locus=", $9); \
            $9 = $9"\;biotype=protein_coding"; \
        } \
        if ($3 == "CDS"){ \
            gsub("product=", "locus=", $9); \
            if ($1 == "MtDNA") gsub("ID=", "ID=CDS:", $9); \
        } \
        if ($3 == "gene"){ \
            $9 = $9"\;biotype=protein_coding"; \
            if ($1 == "MtDNA") gsub("ID=", "ID=gene:", $9); \
        } \
        gsub("hypothetical protein", "hypothetical_protein", $9); \
    }' ${gff} > ${meta.species}.csq.gff
    """

    stub:
    """
    touch ${meta.species}.csf.gff
    """

}