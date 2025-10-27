process LOCAL_VARIANT_MERGE {

    input:
    tuple val(meta), path("annovar.tsv")
    tuple val(meta1), path("csq.tsv")
    tuple val(meta2), path("vep.tsv")
    tuple val(meta3), path("snpeff.tsv")

    output:
    path "variant_transcripts.tsv"

    script:
    """
    cut -f 1,2,10 annovar.tsv | grep -v "N/A" > variant_transcripts.tmp
    cut -f 1,2,10 csq.tsv | grep -v "N/A" >> variant_transcripts.tmp
    cut -f 1,2,10 vep.tsv | grep -v "N/A" >> variant_transcripts.tmp
    cut -f 1,2,9 snpeff.tsv | grep -v "N/A" >> variant_transcripts.tmp
    
    sort -k1,1 -k2,2n -k3,3 variant_transcripts.tmp | uniq > variant_transcripts.tsv

    rm variant_transcripts.tmp    
    """

    stub:
    """
    touch variant_transcripts.tsv
    """

}