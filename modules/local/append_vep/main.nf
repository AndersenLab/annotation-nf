process LOCAL_APPEND_VEP {

    input:
    tuple val(meta), path(tsv)
    tuple val(meta1), path(transcript_tsv)

    output:
    tuple val(meta), path("${meta.species}.${meta.id}.csv.gz")

    script:
    """
    awk 'BEGIN{
        FS="\\t";
        OFS=",";

        grantham["A/A"]="N/A,4";  grantham["A/R"]="112,-1"; grantham["A/N"]="111,-2"; grantham["A/D"]="126,-2"; grantham["A/C"]="195,0";
        grantham["A/Q"]="91,-1";  grantham["A/E"]="107,-1"; grantham["A/G"]="60,0";   grantham["A/H"]="86,-2";  grantham["A/I"]="94,-1";
        grantham["A/L"]="96,-1";  grantham["A/K"]="106,-1"; grantham["A/M"]="84,-1";  grantham["A/F"]="113,-2"; grantham["A/P"]="27,-1";
        grantham["A/S"]="99,1";   grantham["A/T"]="58,0";   grantham["A/W"]="148,-3"; grantham["A/Y"]="112,-2"; grantham["A/V"]="64,0";

        grantham["R/A"]="112,-1"; grantham["R/R"]="N/A,5";  grantham["R/N"]="86,0";   grantham["R/D"]="96,-2";  grantham["R/C"]="180,-3";
        grantham["R/Q"]="43,1";   grantham["R/E"]="54,0";   grantham["R/G"]="125,-2"; grantham["R/H"]="29,0";   grantham["R/I"]="97,-3";
        grantham["R/L"]="102,-2"; grantham["R/K"]="26,2";   grantham["R/M"]="91,-1";  grantham["R/F"]="97,-3";  grantham["R/P"]="103,-2";
        grantham["R/S"]="110,-1"; grantham["R/T"]="71,-1";  grantham["R/W"]="101,-3"; grantham["R/Y"]="77,-2";  grantham["R/V"]="96,-3";

        grantham["N/A"]="111,-2"; grantham["N/R"]="86,0";   grantham["N/N"]="N/A,6";  grantham["N/D"]="23,1";   grantham["N/C"]="139,-3";
        grantham["N/Q"]="46,0";   grantham["N/E"]="42,0";   grantham["N/G"]="80,0";   grantham["N/H"]="68,1";   grantham["N/I"]="149,-3";
        grantham["N/L"]="153,-3"; grantham["N/K"]="94,0";   grantham["N/M"]="142,-2"; grantham["N/F"]="158,-3"; grantham["N/P"]="91,-2";
        grantham["N/S"]="46,1";   grantham["N/T"]="65,0";   grantham["N/W"]="174,-4"; grantham["N/Y"]="143,-2"; grantham["N/V"]="133,-3";

        grantham["D/A"]="126,-2"; grantham["D/R"]="96,-2";  grantham["D/N"]="23,1";   grantham["D/D"]="N/A,6";  grantham["D/C"]="154,-3";
        grantham["D/Q"]="61,0";   grantham["D/E"]="45,2";   grantham["D/G"]="94,-1";  grantham["D/H"]="81,-1";  grantham["D/I"]="168,-3";
        grantham["D/L"]="172,-4"; grantham["D/K"]="101,-1"; grantham["D/M"]="160,-3"; grantham["D/F"]="177,-3"; grantham["D/P"]="108,-1";
        grantham["D/S"]="65,0";   grantham["D/T"]="85,-1";  grantham["D/W"]="181,4";  grantham["D/Y"]="160,-3"; grantham["D/V"]="152,-3";

        grantham["C/A"]="195,0";  grantham["C/R"]="180,-3"; grantham["C/N"]="139,-3"; grantham["C/D"]="154,-3"; grantham["C/C"]="N/A,9";
        grantham["C/Q"]="154,-3"; grantham["C/E"]="158,-4"; grantham["C/G"]="159,-3"; grantham["C/H"]="174,-3"; grantham["C/I"]="198,-1";
        grantham["C/L"]="198,-1"; grantham["C/K"]="202,-3"; grantham["C/M"]="196,-1"; grantham["C/F"]="205,-2"; grantham["C/P"]="169,-3";
        grantham["C/S"]="112,-1"; grantham["C/T"]="149,-1"; grantham["C/W"]="215,-2"; grantham["C/Y"]="194,-2"; grantham["C/V"]="192,-1";

        grantham["Q/A"]="91,-1";  grantham["Q/R"]="43,1";   grantham["Q/N"]="46,0";   grantham["Q/D"]="61,0";   grantham["Q/C"]="154,-3";
        grantham["Q/Q"]="N/A,9";  grantham["Q/E"]="29,2";   grantham["Q/G"]="87,-2";  grantham["Q/H"]="24,0";   grantham["Q/I"]="109,-3";
        grantham["Q/L"]="113,-2"; grantham["Q/K"]="53,1";   grantham["Q/M"]="101,0";  grantham["Q/F"]="116,-3"; grantham["Q/P"]="76,-1";
        grantham["Q/S"]="68,0";   grantham["Q/T"]="42,-1";  grantham["Q/W"]="130,-2"; grantham["Q/Y"]="99,-1";  grantham["Q/V"]="96,-2";

        grantham["E/A"]="107,-1"; grantham["E/R"]="54,0";   grantham["E/N"]="42,0";   grantham["E/D"]="45,2";   grantham["E/C"]="158,-4";
        grantham["E/Q"]="29,2";   grantham["E/E"]="N/A,5";  grantham["E/G"]="98,-2";  grantham["E/H"]="40,0";   grantham["E/I"]="134,-3";
        grantham["E/L"]="138,-3"; grantham["E/K"]="56,1";   grantham["E/M"]="126,-2"; grantham["E/F"]="140,-3"; grantham["E/P"]="93,-1";
        grantham["E/S"]="80,0";   grantham["E/T"]="65,-1";  grantham["E/W"]="152,-3"; grantham["E/Y"]="122,-2"; grantham["E/V"]="121,-2";

        grantham["G/A"]="60,0";   grantham["G/R"]="125,-2"; grantham["G/N"]="80,0";   grantham["G/D"]="94,-1";  grantham["G/C"]="159,-3";
        grantham["G/Q"]="87,-2";  grantham["G/E"]="98,-2";  granthan["G/G"]="N/A,6";  grantham["G/H"]="98,-2";  grantham["G/I"]="135,-4";
        grantham["G/L"]="138,-4"; grantham["G/K"]="127,-2"; grantham["G/M"]="127,-3"; grantham["G/F"]="153,-3"; grantham["G/P"]="42,-2";
        grantham["G/S"]="56,0";   grantham["G/T"]="59,-2";  grantham["G/W"]="184,-2"; grantham["G/Y"]="147,-3"; grantham["G/V"]="109,-3";

        grantham["H/A"]="86,-2";  grantham["H/R"]="29,0";   grantham["H/N"]="68,1";   grantham["H/D"]="81,-1";  grantham["H/C"]="174,-3";
        grantham["H/Q"]="24,0";   grantham["H/E"]="40,0";   grantham["H/G"]="98,-2";  grantham["H/H"]="N/A,8";  grantham["H/I"]="94,-3";
        grantham["H/L"]="99,-3";  grantham["H/K"]="32,-1";  grantham["H/M"]="87,-2";  grantham["H/F"]="100,-1"; grantham["H/P"]="77,-2";
        grantham["H/S"]="89,-1";  grantham["H/T"]="47,-2";  grantham["H/W"]="115,-2"; grantham["H/Y"]="83,-2";  grantham["H/V"]="84,-3";

        grantham["I/A"]="94,-1";  grantham["I/R"]="97,-3";  grantham["I/N"]="149,-3"; grantham["I/D"]="168,-3"; grantham["I/C"]="198,-1";
        grantham["I/Q"]="109,-3"; grantham["I/E"]="134,-3"; grantham["I/G"]="135,-4"; grantham["I/H"]="94,-3";  grantham["I/I"]="N/A,4";
        grantham["I/L"]="5,2";    grantham["I/K"]="102,-3"; grantham["I/M"]="10,1";   grantham["I/F"]="21,0";   grantham["I/P"]="95,-3";
        grantham["I/S"]="142,-2"; grantham["I/T"]="89,-1";  grantham["I/W"]="61,-3";  grantham["I/Y"]="33,-1";  grantham["I/V"]="29,3";

        grantham["L/A"]="96,-1";  grantham["L/R"]="102,-2"; grantham["L/N"]="153,-3"; grantham["L/D"]="172,-4"; grantham["L/C"]="198,-1";
        grantham["L/Q"]="113,-2"; grantham["L/E"]="138,-3"; grantham["L/G"]="138,-4"; grantham["L/H"]="99,-3";  grantham["L/I"]="5,2";
        grantham["L/L"]="N/A,4";  grantham["L/K"]="107,-2"; grantham["L/M"]="15,2";   grantham["L/F"]="22,0";   grantham["L/P"]="98,-3";
        grantham["L/S"]="145,-2"; grantham["L/T"]="92,-1";  grantham["L/W"]="61,-2";  grantham["L/Y"]="36,-1";  grantham["L/V"]="32,1";

        grantham["K/A"]="106,-1"; grantham["K/R"]="26,2";   grantham["K/N"]="94,0";   grantham["K/D"]="101,-1"; grantham["K/C"]="202,-3";
        grantham["K/Q"]="53,1";   grantham["K/E"]="56,1";   grantham["K/G"]="127,-2"; grantham["K/H"]="32,-1";  grantham["K/I"]="102,-3";
        grantham["K/L"]="107,-2"; grantham["K/K"]="N/A,5";  grantham["K/M"]="95,-1";  grantham["K/F"]="102,-3"; grantham["K/P"]="103,-1";
        grantham["K/S"]="121,0";  grantham["K/T"]="78,-1";  grantham["K/W"]="110,-3"; grantham["K/Y"]="85,-2";  grantham["K/V"]="97,-2";

        grantham["M/A"]="84,-1";  grantham["M/R"]="91,-1";  grantham["M/N"]="142,-2"; grantham["M/D"]="160,-3"; grantham["M/C"]="196,-1";
        grantham["M/Q"]="101,0";  grantham["M/E"]="126,-2"; grantham["M/G"]="127,-3"; grantham["M/H"]="87,-2";  grantham["M/I"]="10,1";
        grantham["M/L"]="15,2";   grantham["M/K"]="95,-1";  grantham["M/M"]="N/A,5";  grantham["M/F"]="28,0";   grantham["M/P"]="87,-2";
        grantham["M/S"]="135,-1"; grantham["M/T"]="81,-1";  grantham["M/W"]="67,-1";  grantham["M/Y"]="36,-1";  grantham["M/V"]="21,1";

        grantham["F/A"]="113,-2"; grantham["F/R"]="97,-3";  grantham["F/N"]="158,-3"; grantham["F/D"]="177,-3"; grantham["F/C"]="205,-2";
        grantham["F/Q"]="116,-3"; grantham["F/E"]="140,-3"; grantham["F/G"]="153,-3"; grantham["F/H"]="100,-1"; grantham["F/I"]="21,0";
        grantham["F/L"]="22,0";   grantham["F/K"]="102,-3"; grantham["F/M"]="28,0";   grantham["F/F"]="N/A,6";  grantham["F/P"]="114,-4";
        grantham["F/S"]="155,-2"; grantham["F/T"]="103,-2"; grantham["F/W"]="40,1";   grantham["F/Y"]="22,3";   grantham["F/V"]="50,-1";

        grantham["P/A"]="27,-1";  grantham["P/R"]="103,-2"; grantham["P/N"]="91,-2";  grantham["P/D"]="108,-1"; grantham["P/C"]="169,-3";
        grantham["P/Q"]="76,-1";  grantham["P/E"]="93,-1";  grantham["P/G"]="42,-2";  grantham["P/H"]="77,-2";  grantham["P/I"]="95,-3";
        grantham["P/L"]="98,-3";  grantham["P/K"]="103,-1"; grantham["P/M"]="87,-2";  grantham["P/F"]="114,-4"; grantham["P/P"]="N/A,7";
        grantham["P/S"]="74,-1";  grantham["P/T"]="38,-1";  grantham["P/W"]="147,-4"; grantham["P/Y"]="110,-3"; grantham["P/V"]="68,-2";

        grantham["S/A"]="99,1";   grantham["S/R"]="110,-1"; grantham["S/N"]="46,1";   grantham["S/D"]="65,0";   grantham["S/C"]="112,-1";
        grantham["S/Q"]="68,0";   grantham["S/E"]="80,0";   grantham["S/G"]="56,0";   grantham["S/H"]="89,-1";  grantham["S/I"]="142,-2";
        grantham["S/L"]="145,-2"; grantham["S/K"]="121,0";  grantham["S/M"]="135,-1"; grantham["S/F"]="155,-2"; grantham["S/P"]="74,-1";
        grantham["S/S"]="N/A,4";  grantham["S/T"]="58,1";   grantham["S/W"]="177,-3"; grantham["S/Y"]="144,-2"; grantham["S/V"]="124,-2";

        grantham["T/A"]="58,0";   grantham["T/R"]="71,-1";  grantham["T/N"]="65,0";   grantham["T/D"]="85,-1";  grantham["T/C"]="149,-1";
        grantham["T/Q"]="42,-1";  grantham["T/E"]="65,-1";  grantham["T/G"]="59,-2";  grantham["T/H"]="47,-2";  grantham["T/I"]="89,-1";
        grantham["T/L"]="92,-1";  grantham["T/K"]="78,-1";  grantham["T/M"]="81,-1";  grantham["T/F"]="103,-2"; grantham["T/P"]="38,-1";
        grantham["T/S"]="58,1";   grantham["T/T"]="N/A,5";  grantham["T/W"]="128,-2"; grantham["T/Y"]="92,-2";  grantham["T/V"]="69,0";

        grantham["W/A"]="148,-3"; grantham["W/R"]="101,-3"; grantham["W/N"]="174,-4"; grantham["W/D"]="181,4";  grantham["W/C"]="215,-2";
        grantham["W/Q"]="130,-2"; grantham["W/E"]="152,-3"; grantham["W/G"]="184,-2"; grantham["W/H"]="115,-2"; grantham["W/I"]="61,-3";
        grantham["W/L"]="61,-2";  grantham["W/K"]="110,-3"; grantham["W/M"]="67,-1";  grantham["W/F"]="40,1";   grantham["W/P"]="147,-4";
        grantham["W/S"]="177,-3"; grantham["W/T"]="128,-2"; grantham["W/W"]="N/A,11"; grantham["W/Y"]="37,2";   grantham["W/V"]="88,-3";

        grantham["Y/A"]="112,-2"; grantham["Y/R"]="77,-2";  grantham["Y/N"]="143,-2"; grantham["Y/D"]="160,-3"; grantham["Y/C"]="194,-2";
        grantham["Y/Q"]="99,-1";  grantham["Y/E"]="122,-2"; grantham["Y/G"]="147,-3"; grantham["Y/H"]="83,-2";  grantham["Y/I"]="33,-1";
        grantham["Y/L"]="36,-1";  grantham["Y/K"]="85,-2";  grantham["Y/M"]="36,-1";  grantham["Y/F"]="22,3";   grantham["Y/P"]="110,-3";
        grantham["Y/S"]="144,-2"; grantham["Y/T"]="92,-2";  grantham["Y/W"]="37,2";   grantham["Y/Y"]="N/A,7";  grantham["Y/V"]="55,-1";

        grantham["V/A"]="64,0";   grantham["V/R"]="96,-3";  grantham["V/N"]="133,-3";  grantham["V/D"]="152,-3"; grantham["V/C"]="192,-1";
        grantham["V/Q"]="96,-2";  grantham["V/E"]="121,-2"; grantham["V/G"]="109,-3";  grantham["V/H"]="84,-3";  grantham["V/I"]="29,3";
        grantham["V/L"]="32,1";   grantham["V/K"]="97,-2";  grantham["V/M"]="21,1";    grantham["V/F"]="50,-1";  grantham["V/P"]="68,-2";
        grantham["V/S"]="124,-2"; grantham["V/T"]="69,0";   grantham["V/W"]="88,-3";   grantham["V/Y"]="55,-1";  grantham["V/V"]="N/A,4";
    }{
        if (NR == FNR){
            PP_MAP[\$1] = \$2;
        } else if (FNR == 1) {
            print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, "GRANTHAM,BLOSUM,percent_protein";
        } else {
            if ( \$5 ~ /missense|frameshift|stop|deletion|synonymous|insertion/ ){
                if (\$5 ~ /missense_variant/) {
                    # Extract amino acid substitution from column 7: e.g.  D/N
                    key = \$7;
                    GRANTHAM = grantham[key];
                } else {
                    GRANTHAM = "N/A,N/A";
                }
                key = \$1 "_" \$2 "_" \$10;
                if (key in PP_MAP) PP = PP_MAP[key];
                else PP = "N/A";
            } else {
                GRANTHAM = "N/A,N/A";
                PP = "N/A";
            }
            print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, GRANTHAM, PP;
        }
    }' ${transcript_tsv} ${tsv} | gzip -c > ${meta.species}.${meta.id}.csv.gz
    """

    stub:
    """
    touch ${meta.species}.${meta.id}.csv.gz
    """

}