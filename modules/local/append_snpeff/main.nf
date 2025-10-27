process LOCAL_APPEND_SNPEFF {

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

        # 3-letter to 1-letter amino acid codes
        aa3to1["Ala"] = "A"; aa3to1["Arg"] = "R"; aa3to1["Asn"] = "N"; aa3to1["Asp"] = "D";
        aa3to1["Cys"] = "C"; aa3to1["Gln"] = "Q"; aa3to1["Glu"] = "E"; aa3to1["Gly"] = "G";
        aa3to1["His"] = "H"; aa3to1["Ile"] = "I"; aa3to1["Leu"] = "L"; aa3to1["Lys"] = "K";
        aa3to1["Met"] = "M"; aa3to1["Phe"] = "F"; aa3to1["Pro"] = "P"; aa3to1["Ser"] = "S";
        aa3to1["Thr"] = "T"; aa3to1["Trp"] = "W"; aa3to1["Tyr"] = "Y"; aa3to1["Val"] = "V";

        grantham["A/A"]="N/A"; grantham["A/R"]="112"; grantham["A/N"]="111"; grantham["A/D"]="126"; grantham["A/C"]="195";
        grantham["A/Q"]="91";  grantham["A/E"]="107"; grantham["A/G"]="60";  grantham["A/H"]="86";  grantham["A/I"]="94";
        grantham["A/L"]="96";  grantham["A/K"]="106"; grantham["A/M"]="84";  grantham["A/F"]="113"; grantham["A/P"]="27";
        grantham["A/S"]="99";  grantham["A/T"]="58";  grantham["A/W"]="148"; grantham["A/Y"]="112"; grantham["A/V"]="64";

        grantham["R/A"]="112"; grantham["R/R"]="N/A"; grantham["R/N"]="86";  grantham["R/D"]="96";  grantham["R/C"]="180";
        grantham["R/Q"]="43";  grantham["R/E"]="54";  grantham["R/G"]="125"; grantham["R/H"]="29";  grantham["R/I"]="97";
        grantham["R/L"]="102"; grantham["R/K"]="26";  grantham["R/M"]="91";  grantham["R/F"]="97";  grantham["R/P"]="103";
        grantham["R/S"]="110"; grantham["R/T"]="71";  grantham["R/W"]="101"; grantham["R/Y"]="77";  grantham["R/V"]="96";

        grantham["N/A"]="111"; grantham["N/R"]="86";  grantham["N/N"]="N/A"; grantham["N/D"]="23";  grantham["N/C"]="139";
        grantham["N/Q"]="46";  grantham["N/E"]="42";  grantham["N/G"]="80";  grantham["N/H"]="68";  grantham["N/I"]="149";
        grantham["N/L"]="153"; grantham["N/K"]="94";  grantham["N/M"]="142"; grantham["N/F"]="158"; grantham["N/P"]="91";
        grantham["N/S"]="46";  grantham["N/T"]="65";  grantham["N/W"]="174"; grantham["N/Y"]="143"; grantham["N/V"]="133";

        grantham["D/A"]="126"; grantham["D/R"]="96";  grantham["D/N"]="23";  grantham["D/D"]="N/A"; grantham["D/C"]="154";
        grantham["D/Q"]="61";  grantham["D/E"]="45";  grantham["D/G"]="94";  grantham["D/H"]="81";  grantham["D/I"]="168";
        grantham["D/L"]="172"; grantham["D/K"]="101"; grantham["D/M"]="160"; grantham["D/F"]="177"; grantham["D/P"]="108";
        grantham["D/S"]="65";  grantham["D/T"]="85";  grantham["D/W"]="181"; grantham["D/Y"]="160"; grantham["D/V"]="152";

        grantham["C/A"]="195"; grantham["C/R"]="180"; grantham["C/N"]="139"; grantham["C/D"]="154"; grantham["C/C"]="N/A";
        grantham["C/Q"]="154"; grantham["C/E"]="158"; grantham["C/G"]="159"; grantham["C/H"]="174"; grantham["C/I"]="198";
        grantham["C/L"]="198"; grantham["C/K"]="202"; grantham["C/M"]="196"; grantham["C/F"]="205"; grantham["C/P"]="169";
        grantham["C/S"]="112"; grantham["C/T"]="149"; grantham["C/W"]="215"; grantham["C/Y"]="194"; grantham["C/V"]="192";

        grantham["Q/A"]="91";  grantham["Q/R"]="43";  grantham["Q/N"]="46";  grantham["Q/D"]="61";  grantham["Q/C"]="154";
        grantham["Q/Q"]="N/A"; grantham["Q/E"]="29";  grantham["Q/G"]="87";  grantham["Q/H"]="24";  grantham["Q/I"]="109";
        grantham["Q/L"]="113"; grantham["Q/K"]="53";  grantham["Q/M"]="101"; grantham["Q/F"]="116"; grantham["Q/P"]="76";
        grantham["Q/S"]="68";  grantham["Q/T"]="42";  grantham["Q/W"]="130"; grantham["Q/Y"]="99";  grantham["Q/V"]="96";

        grantham["E/A"]="107"; grantham["E/R"]="54";  grantham["E/N"]="42";  grantham["E/D"]="45";  grantham["E/C"]="158";
        grantham["E/Q"]="29";  grantham["E/E"]="N/A"; grantham["E/G"]="98";  grantham["E/H"]="40";  grantham["E/I"]="134";
        grantham["E/L"]="138"; grantham["E/K"]="56";  grantham["E/M"]="126"; grantham["E/F"]="140"; grantham["E/P"]="93";
        grantham["E/S"]="80";  grantham["E/T"]="65";  grantham["E/W"]="152"; grantham["E/Y"]="122"; grantham["E/V"]="121";

        grantham["G/A"]="60";  grantham["G/R"]="125"; grantham["G/N"]="80";  grantham["G/D"]="94";  grantham["G/C"]="159";
        grantham["G/Q"]="87";  grantham["G/E"]="98";  granthan["G/G"]="N/A"; grantham["G/H"]="98";  grantham["G/I"]="135";
        grantham["G/L"]="138"; grantham["G/K"]="127"; grantham["G/M"]="127"; grantham["G/F"]="153"; grantham["G/P"]="42";
        grantham["G/S"]="56";  grantham["G/T"]="59";  grantham["G/W"]="184"; grantham["G/Y"]="147"; grantham["G/V"]="109";

        grantham["H/A"]="86";  grantham["H/R"]="29";  grantham["H/N"]="68";  grantham["H/D"]="81";  grantham["H/C"]="174";
        grantham["H/Q"]="24";  grantham["H/E"]="40";  grantham["H/G"]="98";  grantham["H/H"]="N/A"; grantham["H/I"]="94";
        grantham["H/L"]="99";  grantham["H/K"]="32";  grantham["H/M"]="87";  grantham["H/F"]="100"; grantham["H/P"]="77";
        grantham["H/S"]="89";  grantham["H/T"]="47";  grantham["H/W"]="115"; grantham["H/Y"]="83";  grantham["H/V"]="84";

        grantham["I/A"]="94";  grantham["I/R"]="97";  grantham["I/N"]="149"; grantham["I/D"]="168"; grantham["I/C"]="198";
        grantham["I/Q"]="109"; grantham["I/E"]="134"; grantham["I/G"]="135"; grantham["I/H"]="94";  grantham["I/I"]="N/A";
        grantham["I/L"]="5";   grantham["I/K"]="102"; grantham["I/M"]="10";  grantham["I/F"]="21";  grantham["I/P"]="95";
        grantham["I/S"]="142"; grantham["I/T"]="89";  grantham["I/W"]="61";  grantham["I/Y"]="33";  grantham["I/V"]="29";

        grantham["L/A"]="96";  grantham["L/R"]="102"; grantham["L/N"]="153"; grantham["L/D"]="172"; grantham["L/C"]="198";
        grantham["L/Q"]="113"; grantham["L/E"]="138"; grantham["L/G"]="138"; grantham["L/H"]="99";  grantham["L/I"]="5";
        grantham["L/L"]="N/A"; grantham["L/K"]="107"; grantham["L/M"]="15";  grantham["L/F"]="22";  grantham["L/P"]="98";
        grantham["L/S"]="145"; grantham["L/T"]="92";  grantham["L/W"]="61";  grantham["L/Y"]="36";  grantham["L/V"]="32";

        grantham["K/A"]="106"; grantham["K/R"]="26";  grantham["K/N"]="94";  grantham["K/D"]="101"; grantham["K/C"]="202";
        grantham["K/Q"]="53";  grantham["K/E"]="56";  grantham["K/G"]="127"; grantham["K/H"]="32";  grantham["K/I"]="102";
        grantham["K/L"]="107"; grantham["K/K"]="N/A"; grantham["K/M"]="95";  grantham["K/F"]="102"; grantham["K/P"]="103";
        grantham["K/S"]="121"; grantham["K/T"]="78";  grantham["K/W"]="110"; grantham["K/Y"]="85";  grantham["K/V"]="97";

        grantham["M/A"]="84";  grantham["M/R"]="91";  grantham["M/N"]="142"; grantham["M/D"]="160"; grantham["M/C"]="196";
        grantham["M/Q"]="101"; grantham["M/E"]="126"; grantham["M/G"]="127"; grantham["M/H"]="87";  grantham["M/I"]="10";
        grantham["M/L"]="15";  grantham["M/K"]="95";  grantham["M/M"]="N/A"; grantham["M/F"]="28";  grantham["M/P"]="87";
        grantham["M/S"]="135"; grantham["M/T"]="81";  grantham["M/W"]="67";  grantham["M/Y"]="36";  grantham["M/V"]="21";

        grantham["F/A"]="113"; grantham["F/R"]="97";  grantham["F/N"]="158"; grantham["F/D"]="177"; grantham["F/C"]="205";
        grantham["F/Q"]="116"; grantham["F/E"]="140"; grantham["F/G"]="153"; grantham["F/H"]="100"; grantham["F/I"]="21";
        grantham["F/L"]="22";  grantham["F/K"]="102"; grantham["F/M"]="28";  grantham["F/F"]="N/A"; grantham["F/P"]="114";
        grantham["F/S"]="155"; grantham["F/T"]="103"; grantham["F/W"]="40";  grantham["F/Y"]="22";  grantham["F/V"]="50";

        grantham["P/A"]="27";  grantham["P/R"]="103"; grantham["P/N"]="91";  grantham["P/D"]="108"; grantham["P/C"]="169";
        grantham["P/Q"]="76";  grantham["P/E"]="93";  grantham["P/G"]="42";  grantham["P/H"]="77";  grantham["P/I"]="95";
        grantham["P/L"]="98";  grantham["P/K"]="103"; grantham["P/M"]="87";  grantham["P/F"]="114"; grantham["P/P"]="N/A";
        grantham["P/S"]="74";  grantham["P/T"]="38";  grantham["P/W"]="147"; grantham["P/Y"]="110"; grantham["P/V"]="68";

        grantham["S/A"]="99";  grantham["S/R"]="110"; grantham["S/N"]="46";  grantham["S/D"]="65";  grantham["S/C"]="112";
        grantham["S/Q"]="68";  grantham["S/E"]="80";  grantham["S/G"]="56";  grantham["S/H"]="89";  grantham["S/I"]="142";
        grantham["S/L"]="145"; grantham["S/K"]="121"; grantham["S/M"]="135"; grantham["S/F"]="155"; grantham["S/P"]="74";
        grantham["S/S"]="N/A"; grantham["S/T"]="58";  grantham["S/W"]="177"; grantham["S/Y"]="144"; grantham["S/V"]="124";

        grantham["T/A"]="58";  grantham["T/R"]="71";  grantham["T/N"]="65";  grantham["T/D"]="85";  grantham["T/C"]="149";
        grantham["T/Q"]="42";  grantham["T/E"]="65";  grantham["T/G"]="59";  grantham["T/H"]="47";  grantham["T/I"]="89";
        grantham["T/L"]="92";  grantham["T/K"]="78";  grantham["T/M"]="81";  grantham["T/F"]="103"; grantham["T/P"]="38";
        grantham["T/S"]="58";  grantham["T/T"]="N/A"; grantham["T/W"]="128"; grantham["T/Y"]="92";  grantham["T/V"]="69";

        grantham["W/A"]="148"; grantham["W/R"]="101"; grantham["W/N"]="174"; grantham["W/D"]="181"; grantham["W/C"]="215";
        grantham["W/Q"]="130"; grantham["W/E"]="152"; grantham["W/G"]="184"; grantham["W/H"]="115"; grantham["W/I"]="61";
        grantham["W/L"]="61";  grantham["W/K"]="110"; grantham["W/M"]="67";  grantham["W/F"]="40";  grantham["W/P"]="147";
        grantham["W/S"]="177"; grantham["W/T"]="128"; grantham["W/W"]="N/A"; grantham["W/Y"]="37";  grantham["W/V"]="88";

        grantham["Y/A"]="112"; grantham["Y/R"]="77";  grantham["Y/N"]="143"; grantham["Y/D"]="160"; grantham["Y/C"]="194";
        grantham["Y/Q"]="99";  grantham["Y/E"]="122"; grantham["Y/G"]="147"; grantham["Y/H"]="83";  grantham["Y/I"]="33";
        grantham["Y/L"]="36";  grantham["Y/K"]="85";  grantham["Y/M"]="36";  grantham["Y/F"]="22";  grantham["Y/P"]="110";
        grantham["Y/S"]="144"; grantham["Y/T"]="92";  grantham["Y/W"]="37";  grantham["Y/Y"]="N/A"; grantham["Y/V"]="55";

        grantham["V/A"]="64";  grantham["V/R"]="96";  grantham["V/N"]="133"; grantham["V/D"]="152"; grantham["V/C"]="192";
        grantham["V/Q"]="96";  grantham["V/E"]="121"; grantham["V/G"]="109"; grantham["V/H"]="84";  grantham["V/I"]="29";
        grantham["V/L"]="32";  grantham["V/K"]="97";  grantham["V/M"]="21";  grantham["V/F"]="50";  grantham["V/P"]="68";
        grantham["V/S"]="124"; grantham["V/T"]="69";  grantham["V/W"]="88";  grantham["V/Y"]="55";  grantham["V/V"]="N/A";
    }{
        if (NR == FNR){
            PP_MAP[\$1] = \$2;
        } else if (FNR == 1) {
            print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, "GRANTHAM,percent_protein";
        } else {
            if ( \$5 ~ /missense|frameshift|stop|deletion|synonymous|insertion/ ){
                if (\$5 ~ /missense/) {
                    # Extract amino acid substitution from column 7: e.g.  p.Ile22Val
                    match(\$7, /^p\\.[A-Z][a-z]{2}/);
                    ref = aa3to1[substr(\$7, RSTART + 2, RLENGTH - 2)];
                    match(\$7, /[A-Z][a-z]{2}\$/);
                    alt = aa3to1[substr(\$7, RSTART, RLENGTH)];
                    key = ref "/" alt;
                    GRANTHAM = grantham[key];
                } else {
                    GRANTHAM = "N/A";
                }
                key = \$1 "_" \$2 "_" \$9;
                if (key in PP_MAP) PP = PP_MAP[key];
                else PP = "N/A";
            } else {
                GRANTHAM = "N/A";
                PP = "N/A";
            }
            print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, GRANTHAM, PP;
        }
    }' ${transcript_tsv} ${tsv} | gzip -c > ${meta.species}.${meta.id}.csv.gz
    """

    stub:
    """
    touch ${meta.species}.${meta.id}.csv.gz
    """

}