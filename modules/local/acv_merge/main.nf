process LOCAL_ACV_MERGE {

    input:
    tuple val(meta), path(tsv1), path(tsv2)
    tuple val(meta1), path(gff_tsv)

    output:
    tuple val(meta), path("${meta.species}.${meta.id}.tsv")

    script:
    """
    #(cat $tsv1; tail -n +2 $tsv2) | sort -t\$'\\t' -k10,10 > sorted.tsv

    #join -t\$'\\t' -1 1 -2 10 -a 2 ${gff_tsv} sorted.tsv | \\
    #    awk -F'\\t' -v OFS='\\t' '{
    #        if (NF == 12) {
    #            print \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$1, \$2, \$3
    #        } else {
    #            print \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$1, \$1, \$1
    #        }
    #    }' > joined.tsv

    #sort -k1,1 -k2,2n joined.tsv > ${meta.species}.${meta.id}.tsv
    #rm sorted.tsv joined.tsv

    (cat $tsv1; tail -n +2 $tsv2) | sort -k1,1 -k2,2n > sorted.tsv

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
                    if (\$10 in transcripts) {
                        TRANSCRIPT = \$10;
                        split(transcripts[\$10], fields, " ");
                        GENE = fields[1];
                        NAME = fields[2];
                    } else if (\$10 in genes) {
                        TRANSCRIPT = "N/A";
                        GENE = \$10;
                        NAME = genes[\$10];
                    } else {
                        TRANSCRIPT = "N/A";
                        GENE = "N/A";
                        NAME = "N/A";
                    }
                }
                print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, TRANSCRIPT, GENE, NAME;
            }
        }' ${gff_tsv} sorted.tsv > ${meta.species}.${meta.id}.tsv

    rm sorted.tsv
    """

    stub:
    """
    touch ${meta.species}.${meta.id}.tsv
    """

}