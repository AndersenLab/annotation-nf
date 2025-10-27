process LOCAL_SNPEFF_MERGE {

    input:
    tuple val(meta), path("raw.tsv")
    tuple val(meta1), path(gff_tsv)

    output:
    tuple val(meta), path("${meta.species}.${meta.id}.tsv")

    script:
    """
    sort -k1,1 -k2,2n raw.tsv > sorted.tsv

    awk -F'\\t' -v OFS='\\t' '{
            if (NR == FNR) {
                genes[\$2] = \$3;
                transcripts[\$1] = \$2 " " \$3;
            } else {
                if (FNR == 1) {
                    TRANSCRIPT = "TRANSCRIPT_NAME";
                    GENE = "WBGENE";
                    NAME = "GENE_NAME";
                } else {
                    if (\$9 in transcripts) {
                        TRANSCRIPT = \$9;
                        split(transcripts[\$9], fields, " ");
                        GENE = fields[1];
                        NAME = fields[2];
                    } else if (\$9 in genes) {
                        TRANSCRIPT = "N/A";
                        GENE = \$9;
                        NAME = genes[\$9];
                    } else {
                        TRANSCRIPT = "N/A";
                        GENE = "N/A";
                        NAME = "N/A";
                    }
                }
                print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, TRANSCRIPT, GENE, NAME;
            }
        }' ${gff_tsv} sorted.tsv > ${meta.species}.${meta.id}.tsv

    rm sorted.tsv
    """

    stub:
    """
    touch ${meta.species}.${meta.id}.tsv
    """

}