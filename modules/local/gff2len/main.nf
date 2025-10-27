process LOCAL_GFF2LEN {

    input:
    tuple val(meta), path(gff)
    path variants

    output:
    tuple val(meta), path("lengths.tsv")

    script:
    """
    awk '
        function binsearch(starts, count, value){
            start = 1;
            end = count + 1;
            while (start < end - 1) {
                mid = int((start + end) / 2);
                if (value >= starts[mid]) {
                    start = mid;
                } else {
                    end = mid;
                }
            }
            return start;
        }
        
        BEGIN{
            FS="\\t";
            OFS="\\t";
        }{

            # If this is the first file, load exon data
            if (FNR == NR){
                if (\$3 == "CDS"){

                    # Get the transcript name
                    match(\$9, /Parent=transcript:[^;]*/);
                    transcript = substr(\$9, RSTART + 18, RLENGTH - 18);  # Skip "Parent=transcript:"

                    # If this is the first time seeing the transcript, initialize array values
                    if (! (transcript in TRANSCRIPT_START)){
                        TRANSCRIPT_START[transcript] = \$4;
                        TRANSCRIPT_STOP[transcript] = \$5;
                        STARTS[transcript] = \$4;
                        STOPS[transcript] = \$5;
                        CUMLENGTH[transcript] = "0";
                        TOTALLENGTH[transcript] = \$5 - \$4 + 1;  # Include both endpoints
                        COUNT[transcript] = 1;

                    # If we have seen this transcript before, add new exon data
                    } else {
                        TRANSCRIPT_STOP[transcript] = \$5;
                        STARTS[transcript] = STARTS[transcript] "," \$4;
                        STOPS[transcript] = STOPS[transcript] "," \$5;
                        CUMLENGTH[transcript] = CUMLENGTH[transcript] "," TOTALLENGTH[transcript];
                        TOTALLENGTH[transcript] = TOTALLENGTH[transcript] + \$5 - \$4 + 1;
                        COUNT[transcript] = COUNT[transcript] + 1;
                    }
                }

            # Otherwise this is the variant file
            } else {

                # Check if variant falls inside associated transcript boundary and proceed if it does
                if (\$3 in TRANSCRIPT_START && \$2 >= TRANSCRIPT_START[\$3] && \$2 <= TRANSCRIPT_STOP[\$3]){
                    KEY = \$1 "_" \$2 "_" \$3;
                    
                    # If this is a single exon transcript, no need to find which exon the variant falls in
                    if (COUNT[\$3] == 1){
                        PP = (\$2 - TRANSCRIPT_START[\$3]) / TOTALLENGTH[\$3];
                        printf "%s\\t%0.2f\\n", KEY, PP*100;

                    # We need to use the binary search to locate which exon the variant occurs in
                    } else {
                        split(STARTS[\$3], starts, ",");
                        split(STOPS[\$3], stops, ",");
                        split(CUMLENGTH[\$3], cumlength, ",");
                       
                        # Use binary search
                        exon_num = binsearch(starts, COUNT[\$3], \$2);
                       
                        # Check if variant falls within the identified exon
                        if (\$2 >= starts[exon_num] && \$2 <= stops[exon_num]){
                            PP = (\$2 - starts[exon_num] + cumlength[exon_num]) / TOTALLENGTH[\$3];
                            printf "%s\\t%0.2f\\n", KEY, PP*100;
                        }
                    }
                }
            }
        }' ${gff} ${variants} > lengths.tsv
    """

    stub:
    """
    touch lengths.tsv
    """
}