#!/usr/bin/env nextflow 
/*
    Authors:
    - Dan Lu <dan.lu@northwestern.edu>
    - Katie Evans <kathrynevans2015@u.northwestern.edu>
    - Ryan McKeown <ryanmckeown2021@u.northwestern.edu>
*/

nextflow.preview.dsl=2
// NXF_VER=20.01.0" Require later version of nextflow
//assert System.getenv("NXF_VER") == "20.01.0"

params.date = new Date().format( 'yyyyMMdd' )
contigs = Channel.from("I","II","III","IV","V","X")
params.ncsq_param = 224
params.help = null

if (params.debug) {
    params.vcf = "${workflow.projectDir}/test_data/WI.20201230.hard-filter.isotype.vcf.gz"
    //params.sample_sheet = "${workflow.projectDir}/test_data/sample_sheet.tsv"
    params.output = "annotation-${params.date}-debug"
    params.divergent_regions = "${workflow.projectDir}/test_data/divergent_regions_strain.bed"

} else {
    // Read input
    params.vcf = null
    //params.sample_sheet = ""
    params.divergent_regions = null

    // folder for the bam files. currently need to put all bam in the same folder
    params.bam_folder = "/projects/b1059/data/${params.species}/WI/alignments/"

    params.output = "annotation-${params.date}"
}

// set default project and ws build for species
if(params.species == "c_elegans") {
    params.project="PRJNA13758"
    params.ws_build="WS283"
} else if(params.species == "c_briggsae") {
    params.project="QX1410_nanopore"
    params.ws_build="Feb2020"
} else if(params.species == "c_tropicalis") {
    params.project="NIC58_nanopore"
    params.ws_build="June2021"
}

    
/* Keep these static ~ They should not need to be updated */
reference_base = "/projects/b1059/data/"
reference_dir = "${reference_base}/${params.species}/genomes/${params.project}/${params.ws_build}"
params.reference = "${reference_dir}/${params.species}.${params.project}.${params.ws_build}.genome.fa.gz"

params.snpeff_reference = "${params.species}.${params.project}.${params.ws_build}"
params.snpeff_dir = "${reference_dir}/snpeff"
snpeff_config = "${reference_dir}/snpeff/snpEff.config"

params.csq_gff = "${reference_dir}/csq/${params.species}.${params.project}.${params.ws_build}.csq.gff3.gz" 
params.AA_score = "${reference_dir}/csq/${params.species}.${params.project}.${params.ws_build}.AA_Scores.tsv"
params.AA_length = "${reference_dir}/csq/${params.species}.${params.project}.${params.ws_build}.AA_Length.tsv"

if(params.species == "c_elegans") {
    params.dust_bed = "${reference_dir}/lcr/${params.species}.${params.project}.${params.ws_build}.dust.bed.gz"
    params.repeat_masker_bed = "${reference_dir}/lcr/${params.species}.${params.project}.${params.ws_build}.repeat_masker.bed.gz"
    params.gene_names="${reference_dir}/csq/wormbase_name_key.txt"

} else if(params.species == "c_briggsae") {
// THESE ARE NOT THE CORRECT FILES - THIS WAS HOW CB WAS RUN IN THE PAST - NEED TO UPDATE
    params.dust_bed = "/projects/b1059/data/c_tropicalis/WI/divergent_regions/20210901/divergent_regions_strain.bed"
    params.repeat_masker_bed = "/projects/b1059/data/c_tropicalis/WI/divergent_regions/20210901/divergent_regions_strain.bed"

    params.gene_names="${reference_dir}/csq/QX1410.R1.current.geneIDs.txt"
} else if(params.species == "c_tropicalis") {
    params.dust_bed = "/projects/b1059/data/c_tropicalis/WI/divergent_regions/20210901/divergent_regions_strain.bed"
    params.repeat_masker_bed = "/projects/b1059/data/c_tropicalis/WI/divergent_regions/20210901/divergent_regions_strain.bed"
    params.gene_names="${reference_dir}/csq/"
}



// Variant annotation files. The same for debug or normal run. 

// "${params.reference_dir}/csq/${species}.${project}.${ws_build}.csq.gff3.gz" from DEC genomes-nf gives some error for transposons
params.snpeff_vcfanno_config = "${workflow.projectDir}/bin/vcfanno_snpeff.toml"
if(params.species == "c_elegans") {
    params.bcsq_vcfanno_config = "${workflow.projectDir}/bin/vcfanno_bcsq.toml"
} else {
    params.bcsq_vcfanno_config = "${workflow.projectDir}/bin/vcfanno_bcsq_ct_cb.toml"
}

// Note that params.species is set in the config to be c_elegans (default)
if ( params.vcf==null ) error "Parameter --vcf is required. Specify path to the full vcf."
//if ( params.sample_sheet==null ) error "Parameter --sample_sheet is required. It should contain a column of strain names, a column of bam file names and a column of bai file names WITH NO HEADERS. If the bam and bai column do not contain full path to the files, specify that path with --bam_folder."
if ( params.species==null ) error "Parameter --species is required. Please select c_elegans, c_briggsae, or c_tropicalis."
if ( params.divergent_regions==null ) println "Parameter --divergent_regions is ignored."


def log_summary() {
/*
    Generates a log
*/

out = '''

Subset isotype reference strains from hard-filter vcf, build trees, define haplotypes and divergent regions.
Should work for c.e, c.b, c.t since they all have the same number of chromosomes.

''' + """
    
-------------    
ANNOTATION-NF
-------------

nextflow main.nf --debug

nextflow main.nf --vcf=hard-filtered.vcf --species=c_elegans --divergent_regions=divergent_regions_strain.bed

    parameters           description                                              Set/Default
    ==========           ===========                                              ========================
    --debug              Set to 'true' to test                                    ${params.debug}
    --species            Species: 'c_elegans', 'c_tropicalis' or 'c_briggsae'     ${params.species}
    --vcf                hard filtered vcf to calculate variant density           ${params.vcf}
    --divergent_regions  Divergent region bed file                                ${params.divergent_regions}
    --reference          Reference used based on species and project              ${params.reference}
    --output             (Optional) output folder name                            ${params.output}
 
    username                                                                      ${"whoami".execute().in.text}

    HELP: http://andersenlab.org/dry-guide/pipeline-annotation-nf   
    ----------------------------------------------------------------------------------------------
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId] 
"""
out
}


log.info(log_summary())


if (params.help) {
    exit 1
}


// Read input
input_vcf = Channel.fromPath("${params.vcf}")
input_vcf_index = Channel.fromPath("${params.vcf}.tbi")

// To convert ref strain to isotype names. Note only strains that match the first col will get changed, so won't impact ct and cb.
isotype_convert_table = Channel.fromPath("${workflow.projectDir}/bin/ref_strain_isotype.tsv") 

workflow { 

    // snpeff annotation
    input_vcf.combine(input_vcf_index)
      .combine(Channel.fromPath(params.snpeff_vcfanno_config))
      .combine(Channel.fromPath(params.dust_bed))
      .combine(Channel.fromPath(params.dust_bed + ".tbi"))
      .combine(Channel.fromPath(params.repeat_masker_bed))
      .combine(Channel.fromPath(params.repeat_masker_bed + ".tbi")) | snpeff_annotate_vcf

    
    // bcsq annotation
    input_vcf.combine(input_vcf_index)
      .combine(Channel.fromPath(params.csq_gff)) | bcsq_annotate_vcf

    // add blossum and perc protein etc. to vcf
    bcsq_annotate_vcf.out.bcsq_tsv
      .combine(Channel.fromPath(params.AA_length))
      .combine(Channel.fromPath(params.AA_score)) | prep_other_annotation


    // if not celegans, use different vcfanno file without dust bed
    bcsq_annotate_vcf.out.bcsq_vcf
      .combine(prep_other_annotation.out)
      .combine(Channel.fromPath(params.bcsq_vcfanno_config))
      .combine(Channel.fromPath(params.dust_bed))
      .combine(Channel.fromPath(params.dust_bed + ".tbi"))
      .combine(Channel.fromPath(params.repeat_masker_bed))
      .combine(Channel.fromPath(params.repeat_masker_bed + ".tbi")) | AA_annotate_vcf
    

    // parse bcsq annotation to make flat file
    AA_annotate_vcf.out.vcfanno_vcf | bcsq_extract_scores | bcsq_parse_scores

    if(params.divergent_regions) {
        AA_annotate_vcf.out.vcfanno_vcf
            .combine(Channel.fromPath(params.divergent_regions)) | bcsq_extract_samples | bcsq_parse_samples
    } else {
        AA_annotate_vcf.out.vcfanno_vcf | bcsq_extract_samples | bcsq_parse_samples
    }
    

    bcsq_parse_samples.out
    	.combine(bcsq_parse_scores.out)
        .combine(Channel.fromPath(params.gene_names))
        .combine(snpeff_annotate_vcf.out.snpeff_flat) | make_flat_file
    	// .combine(Channel.fromPath("${workflow.projectDir}/bin/wormbase_name_key.txt")) | make_flat_file

      
	// Generate Strain-level TSV and VCFs. 
	// note 1: since annotation is only done for isotype-ref strains, here also only include isotype ref strains.
	// note 2: this step used vcf with only bcsq annotation b/c I'm not sure how to split into single sample vcf with protein length and amino acid score. Ryan has the code to split out single strain annotation for bcsq.
	  
	bcsq_annotate_vcf.out.bcsq_vcf | strain_list
	strain_set = strain_list.out.splitCsv(sep: ' ', strip: true)

	strain_set.combine( AA_annotate_vcf.out.vcfanno_vcf ) | generate_strain_vcf
	  
	// // Extract severity tracks
	snpeff_tracks = Channel.from(["LOW", "MODERATE", "HIGH", "MODIFIER"])  
    bcsq_tracks = Channel.from(["LOW", "HIGH"])

    // dont need these anymore
    /*
	snpeff_annotate_vcf.out.snpeff_vcf.combine(snpeff_tracks) | snpeff_severity_tracks
	make_flat_file.out.combine(bcsq_tracks) | bcsq_severity_tracks
    */

}


/* 
    =========================
    Annotate isotype-only vcf
    =========================
*/


process snpeff_annotate_vcf {

    // conda "/projects/b1059/software/conda_envs/popgen-nf_env"
    label 'annotation'

    publishDir "${params.output}/variation", mode: 'copy', pattern: '*.snpeff.vcf*'
    publishDir "${params.output}/variation", mode: 'copy', pattern: 'snpeff.stats.csv'

    input:
        tuple file(vcf), file(vcf_index), \
              path(vcfanno), \
              path("dust.bed.gz"), \
              path("dust.bed.gz.tbi"), \
              path("repeat_masker.bed.gz"), \
              path("repeat_masker.bed.gz.tbi")

    output:
        tuple path("*.snpeff.vcf.gz"), path("*.snpeff.vcf.gz.tbi"), emit: snpeff_vcf
        path "snpeff.stats.csv"
        path "WI-snpeff.tsv", emit: snpeff_flat


    script:
    """
        bcftools view -O v ${vcf} | \\
        snpEff eff -csvStats snpeff.stats.csv \\
                   -no-downstream \\
                   -no-intergenic \\
                   -no-upstream \\
                   -nodownload \\
        -dataDir ${params.snpeff_dir} \\
        -config ${params.snpeff_dir}/snpEff.config \\
        ${params.snpeff_reference} | \\
        bcftools view -O z > out.vcf.gz


        if [ ${params.species} = "c_elegans" ]
        then
            vcfanno ${vcfanno} out.vcf.gz | bcftools view -O z > WI.${params.date}.hard-filter.isotype.snpeff.vcf.gz
            bcftools index --tbi WI.${params.date}.hard-filter.isotype.snpeff.vcf.gz
        else
            mv out.vcf.gz WI.${params.date}.hard-filter.isotype.snpeff.vcf.gz
            bcftools index --tbi WI.${params.date}.hard-filter.isotype.snpeff.vcf.gz
        fi

        bcftools query -e 'INFO/ANN="."' -f'%CHROM %POS %ID %REF %ALT %QUAL %FILTER %ANN\n' WI.${params.date}.hard-filter.isotype.snpeff.vcf.gz > WI-snpeff.tsv

    """

}



process bcsq_annotate_vcf {

    // bcsq needs bcftools=1.12 to work properly
    // conda "/projects/b1059/software/conda_envs/bcftools"
    label 'annotation'

    memory 16.GB

    input:
        tuple file(vcf), file(vcf_index), file(gff)


    output:
        tuple file("bcsq.vcf.gz"), file("bcsq.vcf.gz.tbi"), emit:bcsq_vcf
        path "WI-BCSQ.tsv", emit:bcsq_tsv

    script:
    """
        # bgzip -c $gff > csq.gff.gz
        cp $gff csq.gff.gz
        # tabix -p gff csq.gff.gz

        bcftools csq -O z --fasta-ref ${params.reference} \\
                     --gff-annot csq.gff.gz \\
                     --ncsq ${params.ncsq_param} \\
                     --phase a $vcf > bcsq.vcf.gz

        bcftools index --tbi bcsq.vcf.gz

        bcftools query -e 'INFO/BCSQ="."' -f'%CHROM %POS %ID %REF %ALT %QUAL %FILTER %BCSQ\n' bcsq.vcf.gz > WI-BCSQ.tsv
    """

}

process prep_other_annotation {

    // conda "/projects/b1059/software/conda_envs/popgen-nf-r_env"
    label 'R'
    echo true
    memory 16.GB

    input:
        tuple path("WI-BCSQ.tsv"), path("gff_AA_Length.tsv"), path("AA_Scores.tsv")

    output:
        path("anno_vcf.vcf")

    """
      Rscript --vanilla ${workflow.projectDir}/bin/make_vcf_annotation.R WI-BCSQ.tsv AA_Scores.tsv gff_AA_Length.tsv
    """

  // check that AA length gets changed with new gff in new wb version (genomes-nf)

}


process AA_annotate_vcf {

    // conda "/projects/b1059/software/conda_envs/popgen-nf_env"
    label 'annotation'

    publishDir "${params.output}/variation", mode: 'copy'

    memory 16.GB

    input:
        tuple file(vcf), file(vcf_index), file("anno_vcf.vcf"), file(vcfanno), \
        file("dust.bed.gz"), file("dust.bed.gz.tbi"), file("repeat_masker.bed.gz"), file("repeat_masker.bed.gz.tbi")

    output:
        tuple file("*hard-filter.isotype.bcsq.vcf.gz"), file("*hard-filter.isotype.bcsq.vcf.gz.tbi"), emit: vcfanno_vcf
        file("*.stats.txt")


    """
        bcftools reheader -h ${workflow.projectDir}/bin/new_header anno_vcf.vcf > header_anno_vcf.vcf
        bgzip header_anno_vcf.vcf
		bcftools index -t header_anno_vcf.vcf.gz

        # change vcfanno config
        cat ${vcfanno} | sed 's/REPLACE_HEADER_VCF/header_anno_vcf.vcf.gz/' > vcfanno.config

		vcfanno vcfanno.config $vcf | bcftools view -O z > WI.${params.date}.hard-filter.isotype.bcsq.vcf.gz
        bcftools index --tbi WI.${params.date}.hard-filter.isotype.bcsq.vcf.gz
        bcftools stats -s- WI.${params.date}.hard-filter.isotype.bcsq.vcf.gz > WI.${params.date}.hard-filter.isotype.bcsq.vcf.gz.stats.txt

    """
}

process bcsq_extract_samples {


    // conda "/projects/b1059/software/conda_envs/bcftools"
    label 'annotation'

    memory 16.GB

    input:
        tuple file(bcsq_vcf), file(bcsq_vcf_index), file(divergent_regions)

    output:
        tuple file("BCSQ_samples.tsv"), file("div.bed")

    """
        # Extract samples
        bcftools view -e 'INFO/BCSQ="."' -Ou ${bcsq_vcf}  | \\
        bcftools query -i 'GT ="alt"' -f'%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE:%TBCSQ{*}=]\n' > BCSQ_samples.tsv

        # format divergent file
        bedtools merge -i ${divergent_regions} > div.bed
    """

}

process bcsq_extract_scores {


    // conda "/projects/b1059/software/conda_envs/bcftools"
    label 'annotation'

    memory 16.GB

    input:
        tuple file(bcsq_vcf), file(bcsq_vcf_index)

    output:
        path "BCSQ_scores.tsv"

    """
        # Extract scores
        bcftools view -e 'INFO/BCSQ="."' -Ou ${bcsq_vcf} | \\
        bcftools query -i 'GT ="alt"' -f'%CHROM\t%POS\t%REF\t%ALT\t%BCSQ\t%BLOSUM\t%GRANTHAM\t%PERCENT_PROTEIN\n' > BCSQ_scores.tsv

    """

}


process bcsq_parse_scores {


    // conda "/projects/b1059/software/conda_envs/popgen-nf-r_env"
    label 'R'

    memory 16.GB

    input:
        file(bcsq_scores)

    output:
        path "BCSQ_scores_parsed.tsv"

    """
        Rscript --vanilla ${workflow.projectDir}/bin/parse_scores_bcsq.R $bcsq_scores
    """

}

process bcsq_parse_samples {


    // conda "/projects/b1059/software/conda_envs/popgen-nf-r_env"
    label 'R'

    memory 64.GB

    input:
        tuple file(bcsq_samples), file(div_file)

    output:
        tuple file("BCSQ_samples_parsed.tsv"), file("div.bed")

    """
        Rscript --vanilla ${workflow.projectDir}/bin/parse_samples_bcsq.R $bcsq_samples
    """

}

process make_flat_file {


    // conda "/projects/b1059/software/conda_envs/popgen-nf-r_env"
    label 'R'
    publishDir "${params.output}/variation", mode: 'copy'

    //memory 16.GB

    memory 48.GB

    input:
        tuple file(bcsq_samples_parsed), file(div_file), file(bcsq_scores), file(wbgene_names), file(snpeff)

    output:
        path "*.strain-annotation.tsv"

    """
        Rscript --vanilla ${workflow.projectDir}/bin/make_flat_file.R $bcsq_samples_parsed $bcsq_scores $wbgene_names $div_file ${snpeff} ${params.species}

        # rename flatfile
        mv WI-annotated-flatfile.tsv WI.${params.date}.strain-annotation.tsv
    """

}


/* 
    ==============================
    Generate individual strain vcf.
    ==============================
*/


process strain_list {

    // conda "/projects/b1059/software/conda_envs/popgen-nf_env"
    label 'annotation'

    input:
        tuple path(vcf), path(vcf_index)
    
    output:
        file("samples.txt")

    """
        bcftools query --list-samples ${vcf} > samples.txt
    """
}


process generate_strain_vcf {
    // Generate a single VCF for every strain.

    // conda "/projects/b1059/software/conda_envs/popgen-nf_env"
    label 'annotation'

    tag { strain }

    publishDir "${params.output}/strain_vcf", mode: 'copy', pattern: "*.vcf.gz*"

    input:
        tuple val(strain), path(vcf), file(vcf_index)

    output:
        path("*.vcf.gz*")
        // tuple path("${strain}.${date}.vcf.gz"),  path("${strain}.${date}.vcf.gz.tbi"), path("${strain}.${date}.vcf.gz.csi")

    """
        bcftools view -s ${strain} ${vcf} > ${strain}.${params.date}.vcf
        bgzip ${strain}.${params.date}.vcf
        bcftools index -c ${strain}.${params.date}.vcf.gz
        bcftools index -t ${strain}.${params.date}.vcf.gz

        # keep both isotype and ref strain
        if grep -q $strain ${workflow.projectDir}/bin/isotype_ref_strain.tsv; 
        then
            new_strain=`grep ${strain} ${workflow.projectDir}/bin/isotype_ref_strain.tsv | awk '{print \$2}'`
            bcftools reheader --samples ${workflow.projectDir}/bin/isotype_ref_strain.tsv ${strain}.${params.date}.vcf.gz | \\
            bcftools view -Oz > \$new_strain.20220207.vcf.gz

            bcftools index -c \$new_strain.20220207.vcf.gz
            bcftools index -t \$new_strain.20220207.vcf.gz
        fi

    """

}


/*
process generate_strain_tsv {
    // Generate a single TSV for every strain.

    tag { strain }

    publishDir "${params.output}/strain/tsv", mode: 'copy', pattern: "*.tsv.gz*"

    input:
        tuple val(strain), path(vcf), file(vcf_index)

    output:
        tuple path("${strain}.${date}.tsv.gz"),  path("${strain}.${date}.tsv.gz.tbi")

    """
        # Generate TSV
        {
            echo -e 'CHROM\\tPOS\\tREF\\tALT\\tFILTER\\tFT\\tGT';
            bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%FT\t%TGT]\n' --samples ${strain} ${vcf};
        } | awk -F'\\t' -vOFS='\\t' '{ gsub("\\\\.", "PASS", \$6) ; print }' > ${strain}.${date}.tsv

        bgzip ${strain}.${date}.tsv
        tabix -S 1 -s 1 -b 2 -e 2 ${strain}.${date}.tsv.gz
    """

}
*/



/* 
    ==========================================
    Generate severity tracks for cendr browser
    ==========================================
*/


process snpeff_severity_tracks {
    /*
        The severity tracks are bedfiles with annotations for
        LOW
        MODERATE
        HIGH
        MODIFIER
        variants as annotated with SNPEff
        They are used on the CeNDR browser. Previous CeNDR browser used soft filtered vcf, this one uses hard filtered vcf
    */

    // conda "/projects/b1059/software/conda_envs/popgen-nf_env"
    label 'annotation'

    publishDir "${params.output}/tracks/snpeff/", mode: 'copy'

    tag { severity }

    input:
        tuple path("in.vcf.gz"), path("in.vcf.gz.csi"), val(severity)
    output:
        set file("${params.date}.${severity}.bed.gz"), file("${params.date}.${severity}.bed.gz.tbi")

    """
        bcftools view --apply-filters PASS in.vcf.gz | \
        grep ${severity} | \
        awk '\$0 !~ "^#" { print \$1 "\\t" (\$2 - 1) "\\t" (\$2)  "\\t" \$1 ":" \$2 "\\t0\\t+"  "\\t" \$2 - 1 "\\t" \$2 "\\t0\\t1\\t1\\t0" }' | \\
        bgzip  > ${params.date}.${severity}.bed.gz
        tabix -p bed ${params.date}.${severity}.bed.gz
        fsize=\$(zcat ${params.date}.${severity}.bed.gz | wc -c)
        if [ \${fsize} -lt 2000 ]; then
            exit 1
        fi;
    """
}


process bcsq_severity_tracks {

    // conda "/projects/b1059/software/conda_envs/popgen-nf_env"
    label 'annotation'

    publishDir "${params.output}/tracks/bcsq/", mode: 'copy'

    tag { severity }

    input:
        tuple path("flatfile"), val(severity)
    output:
        tuple file("${params.date}.${severity}.bed.gz"), file("${params.date}.${severity}.bed.gz.tbi")

    """
        cat ${flatfile} | \
        grep ${severity} | \
        awk '\$0 !~ "^#" { print \$1 "\\t" (\$2 - 1) "\\t" (\$2)  "\\t" \$1 ":" \$2 "\\t0\\t+"  "\\t" \$2 - 1 "\\t" \$2 "\\t0\\t1\\t1\\t0" }' | \\
        bgzip  > ${params.date}.${severity}.bed.gz
        tabix -p bed ${params.date}.${severity}.bed.gz
        fsize=\$(zcat ${params.date}.${severity}.bed.gz | wc -c)
        if [ \${fsize} -lt 2000 ]; then
            exit 1
        fi;
    """
}


