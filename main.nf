#!/usr/bin/env nextflow 

// Needed to publish results
nextflow.preview.output = true

if (params.species != null) {
    if (params.vcf != null){
        if (params.vcf.toString().isInteger()){
            release = params.vcf
            vcf = "${params.data_dir}/${params.species}/WI/variation/${release}/vcf/WI.${params.vcf}.hard-filter.isotype.vcf.gz"
        } else {
            release = null
            vcf = params.vcf
        }
    } else {
        println "A vcf file must be specified using --vcf either as a file path or a valid CAENDR release date"
        exit 1
    }
    if (params.hdr == null) {
        if (params.species == "c_elegans") {
            release ?= "20250625"
        } else if (params.species == "c_briggsae") {
            release ?= "20250626"
        } else if (params.species == "c_tropicalis") {
            release ?= "20250627"
        }
        if (release != null) {
            hdr = "${params.data_dir}/${params.species}/WI/divergent_regions/${release}/${release}_${params.species}_divergent_regions_strain.bed.gz"
        } else {
            println "If no HDR file is specified, a valid CAENDR species must be specified"
        }
    } else {
        hdr = params.hdr
    }
    if (params.gff == null || params.genome == null){
        if (params.species == "c_elegans") {
            project = "PRJNA13758"
            build = "WS283"
        } else if (params.species == "c_briggsae") {
            project = "QX1410_nanopore"
            build = "Feb2020"
        } else if (params.species == "c_tropicalis") {
            project = "NIC58_nanopore"
            build = "June2021"
        } else {
            println "The GFF and genome file paths must be specified using --gff and --genome if the species is not a CAENDR species"
            exit 1
        }
        if (params.gff == null) {
            gff = "${params.data_dir}/${params.species}/genomes/${project}/${build}/csq/${params.species}.${project}.${build}.csq.gff3"
        } else {
            gff = params.gff
        }
        if (params.genome == null) {
            genome = "${params.data_dir}/${params.species}/genomes/${project}/${build}/${params.species}.${project}.${build}.genome.fa"
        } else {
            genome = params.genome
        }
    } else {
        gff = params.gff
        genome = params.genome
    }
} else if (params.help == false) {
    if (params.vcf == null || params.gff == null || params.genome == null || params.hdr == null) {
        println "If no species is specified, valid VCF, GFF, genome, and HDR file paths must be passed using --vcf, --gcf, --genome, and --hdr, respectively"
        exit 1
    } else {
        vcf = params.vcf
        gff = params.gff
        genome = params.genome
        hdr = params.hdr
    }
}


def log_summary() {
/*
    Generates a log
*/

out = '''

Create variant annotations using the tools Annovar, CSQ, SnpEff, and VEP

''' + """
    
-------------    
ANNOTATION-NF
-------------

nextflow main.nf --help

nextflow main.nf --vcf=hard-filtered.vcf --species=c_elegans --divergent_regions=divergent_regions_strain.bed

    parameters           description                                              Set/Default
    ==========           ===========                                              ========================
    --species            Species: 'c_elegans', 'c_tropicalis' or 'c_briggsae'     ${params.species}
    --vcf                hard filtered vcf or release date                        ${vcf}
    --gff                Csq-compatible gff file                                  ${gff}
    --hdr                Divergent region strain bed file                         ${hdr}
    --genome            Reference used based on species and project               ${genome}
    -output-dir          (Optional) output folder name                            ${workflow.outputDir}
 
    username                                                                      ${"whoami".execute().in.text}

    HELP: http://andersenlab.org/dry-guide/pipeline-annotation-nf   
    ----------------------------------------------------------------------------------------------
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId] 
"""
out
}


log.info(log_summary())


include { BCFTOOLS_STRAIN_MATRIX                 } from './modules/bcftools/strain_matrix/main.nf'
include { BEDTOOLS_HDR_ANNOTATIONS               } from './modules/bedtools/hdr_annotations/main.nf'
include { BCFTOOLS_HDR_MATRIX                    } from './modules/bcftools/hdr_matrix/main.nf'
include { LOCAL_GFF2TSV                          } from './modules/local/gff2tsv/main.nf'
include { PYTHON_MT_GFF2GTF                      } from './modules/python/mt_gff2gtf/main.nf'
include { PERL_ANNOVAR_ANNOTATION                } from './modules/perl/annovar_annotation/main.nf'
include { BCFTOOLS_REFORMAT_ANNOVAR              } from './modules/bcftools/reformat_annovar/main.nf'
include { LOCAL_ACV_MERGE as LOCAL_ANNOVAR_MERGE } from './modules/local/acv_merge/main.nf'
include { BCFTOOLS_CSQ_ANNOTATION                } from './modules/bcftools/csq_annotation/main.nf'
include { BCFTOOLS_REFORMAT_CSQ                  } from './modules/bcftools/reformat_csq/main.nf'
include { LOCAL_ACV_MERGE as LOCAL_CSQ_MERGE     } from './modules/local/acv_merge/main.nf'
include { VEP_ANNOTATION                         } from './modules/vep/annotation/main.nf'
include { BCFTOOLS_REFORMAT_VEP                  } from './modules/bcftools/reformat_vep/main.nf'
include { LOCAL_ACV_MERGE as LOCAL_VEP_MERGE     } from './modules/local/acv_merge/main.nf'
include { JAVA_SNPEFF_ANNOTATION                 } from './modules/java/snpeff_annotation/main.nf'
include { BCFTOOLS_REFORMAT_SNPEFF               } from './modules/bcftools/reformat_snpeff/main.nf'
include { LOCAL_SNPEFF_MERGE                     } from './modules/local/snpeff_merge/main.nf'
include { BCFTOOLS_COMPRESS_INDEX                } from './modules/bcftools/compress_index/main.nf'
include { LOCAL_VARIANT_MERGE                    } from './modules/local/variant_merge/main.nf'
include { LOCAL_GFF2LEN                          } from './modules/local/gff2len/main.nf'
include { LOCAL_APPEND_ANNOVAR                   } from './modules/local/append_annovar/main.nf'
include { LOCAL_APPEND_CSQ                       } from './modules/local/append_csq/main.nf'
include { LOCAL_APPEND_VEP                       } from './modules/local/append_vep/main.nf'
include { LOCAL_APPEND_SNPEFF                    } from './modules/local/append_snpeff/main.nf'


if (params.help) {
    exit 1
}



workflow {
    main:
    ch_versions = Channel.empty()

    vcf_ch = Channel.fromPath( vcf, checkIfExists: true ).map { it: [[species: params.species], it, "${it}.tbi"] }
    gff_ch = Channel.fromPath( gff, checkIfExists: true ).map { it: [[species: params.species], it] }
    hdr_ch = Channel.fromPath( hdr, checkIfExists: true ).first()
    genome_ch = Channel.fromPath( genome, checkIfExists: true )
    option_ch = Channel.of( 0, 1 )

    BCFTOOLS_STRAIN_MATRIX( vcf_ch )
    ch_versions = ch_versions.mix(BCFTOOLS_STRAIN_MATRIX.out.versions)
    mt_vcf_ch = BCFTOOLS_STRAIN_MATRIX.out.mt
    no_mt_vcf_ch = BCFTOOLS_STRAIN_MATRIX.out.no_mt
    matrix_ch = BCFTOOLS_STRAIN_MATRIX.out.matrix.first()
    sample_ch = BCFTOOLS_STRAIN_MATRIX.out.samples.splitCsv( strip: true ).map { it: it[0] }

    BEDTOOLS_HDR_ANNOTATIONS( matrix_ch,
                              BCFTOOLS_STRAIN_MATRIX.out.no_mt_bed.first(),
                              hdr_ch,
                              sample_ch )
    ch_versions = ch_versions.mix(BEDTOOLS_HDR_ANNOTATIONS.out.versions.first())
    
    BCFTOOLS_HDR_MATRIX( matrix_ch,
                         BEDTOOLS_HDR_ANNOTATIONS.out.geno.collect(),
                         no_mt_vcf_ch )
    ch_versions = ch_versions.mix(BCFTOOLS_HDR_MATRIX.out.versions)
    hdr_vcf_ch = BCFTOOLS_HDR_MATRIX.out.hdr
    
    LOCAL_GFF2TSV( gff_ch )
    gff_tsv_ch = LOCAL_GFF2TSV.out.first()

    PYTHON_MT_GFF2GTF( gff_ch,
                    Channel.fromPath( "${workflow.projectDir}/bin/mt_gff2gtf.py") )
    ch_versions = ch_versions.mix(PYTHON_MT_GFF2GTF.out.versions)
    mt_gtf_ch = PYTHON_MT_GFF2GTF.out.gtf


    ///////////////
    // ANNOVAR
    ///////////////

    PERL_ANNOVAR_ANNOTATION( hdr_vcf_ch,
                             genome_ch,
                             gff_ch )
    annovar_ch = PERL_ANNOVAR_ANNOTATION.out.annovar.map { it: [[species: it[0].species, id: "annovar"], it[1]] }.first()
    ch_versions = ch_versions.mix(PERL_ANNOVAR_ANNOTATION.out.versions)
    
    BCFTOOLS_REFORMAT_ANNOVAR( annovar_ch,
                               option_ch )
    ch_versions = ch_versions.mix(BCFTOOLS_REFORMAT_ANNOVAR.out.versions.first())

    annovar_flat_ch = BCFTOOLS_REFORMAT_ANNOVAR.out.tsv
        .collect()
        .map { it:
            if (it[1].endsWith('.hdr0.annovar.tsv')) {
                [it[0], it[1], it[3]]
            } else {
                [it[0], it[3], it[1]]
            }
        }

    LOCAL_ANNOVAR_MERGE( annovar_flat_ch,
                         gff_tsv_ch )
    merged_annovar_ch = LOCAL_ANNOVAR_MERGE.out
        

    ///////////////
    // CSQ
    ///////////////

    BCFTOOLS_CSQ_ANNOTATION( hdr_vcf_ch,
                             genome_ch,
                             gff_ch )
    csq_ch = BCFTOOLS_CSQ_ANNOTATION.out.csq.map { it: [[species: it[0].species, id: "csq"], it[1], it[2]]}.first()
    ch_versions = ch_versions.mix(BCFTOOLS_CSQ_ANNOTATION.out.versions)

    BCFTOOLS_REFORMAT_CSQ( csq_ch,
                           option_ch )
    ch_versions = ch_versions.mix(BCFTOOLS_REFORMAT_CSQ.out.versions.first())

    csq_flat_ch = BCFTOOLS_REFORMAT_CSQ.out.tsv
        .collect()
        .map { it:
            if (it[1].endsWith('.hdr0.csq.tsv')) {
                [it[0], it[1], it[3]]
            } else {
                [it[0], it[3], it[1]]
            }
        }
    LOCAL_CSQ_MERGE( csq_flat_ch,
                     gff_tsv_ch )
    merged_csq_ch = LOCAL_CSQ_MERGE.out


    ///////////////
    // VEP
    ///////////////

    VEP_ANNOTATION( hdr_vcf_ch,
                    genome_ch,
                    gff_ch )
    vep_ch = VEP_ANNOTATION.out.vep.map { it: [[species: it[0].species, id: "vep"], it[1]]}.first()
    ch_versions = ch_versions.mix(VEP_ANNOTATION.out.versions)

    BCFTOOLS_REFORMAT_VEP( vep_ch,
                           option_ch )
    ch_versions = ch_versions.mix(BCFTOOLS_REFORMAT_VEP.out.versions.first())

    vep_flat_ch = BCFTOOLS_REFORMAT_VEP.out.tsv
        .collect()
        .map { it:
            if (it[1].endsWith('.hdr0.vep.tsv')) {
                [it[0], it[1], it[3]]
            } else {
                [it[0], it[3], it[1]]
            }
        }
    LOCAL_VEP_MERGE( vep_flat_ch,
                     gff_tsv_ch )
    merged_vep_ch = LOCAL_VEP_MERGE.out


    ///////////////
    // SNPEFF
    ///////////////

    JAVA_SNPEFF_ANNOTATION( mt_vcf_ch,
                            genome_ch,
                            mt_gtf_ch,
                            Channel.fromPath("${workflow.projectDir}/bin/base_snpeff.config") )
    ch_versions = ch_versions.mix(JAVA_SNPEFF_ANNOTATION.out.versions)

    BCFTOOLS_REFORMAT_SNPEFF( JAVA_SNPEFF_ANNOTATION.out.snpeff.map { it: [[species: it[0].species, id: "snpeff"], it[1]]} )
    ch_versions = ch_versions.mix(BCFTOOLS_REFORMAT_SNPEFF.out.versions)

    LOCAL_SNPEFF_MERGE( BCFTOOLS_REFORMAT_SNPEFF.out.tsv,
                        gff_tsv_ch )
    merged_snpeff_ch = LOCAL_SNPEFF_MERGE.out


    //////////////////////////
    // Append variant features
    //////////////////////////

    uncompress_ch = PERL_ANNOVAR_ANNOTATION.out.annovar
        .mix(JAVA_SNPEFF_ANNOTATION.out.snpeff)
        .mix(VEP_ANNOTATION.out.vep)

    BCFTOOLS_COMPRESS_INDEX(uncompressed_ch)
    ch_versions = ch_versions.mix(BCFTOOLS_COMPRESS_INDEX.out.versions.first())

    LOCAL_VARIANT_MERGE( merged_annovar_ch,
                         merged_csq_ch,
                         merged_vep_ch,
                         merged_snpeff_ch )

    LOCAL_GFF2LEN( gff_ch,
                   LOCAL_VARIANT_MERGE.out )
    gff_len_ch = LOCAL_GFF2LEN.out.first()

    LOCAL_APPEND_ANNOVAR( merged_annovar_ch,
                          gff_len_ch )

    LOCAL_APPEND_CSQ( merged_csq_ch,
                      gff_len_ch )

    LOCAL_APPEND_VEP( merged_vep_ch,
                      gff_len_ch )

    LOCAL_APPEND_SNPEFF( merged_snpeff_ch,
                         gff_len_ch )

    // Collate and save software versions
        ch_versions
            .collectFile(name: 'workflow_software_versions.txt', sort: true, newLine: true)
            .set { ch_collated_versions }

    publish:
    LOCAL_APPEND_ANNOVAR.out            >> "annotation"
    LOCAL_APPEND_CSQ.out                >> "annotation"
    BCFTOOLS_CSQ_ANNOTATION.out         >> "variation"
    LOCAL_APPEND_VEP.out                >> "annotation"
    LOCAL_APPEND_SNPEFF.out             >> "annotation"
    BCFTOOLS_COMPRESS_INDEX.out.indexed >> "variation"
    ch_collated_versions                >> "."
}

// Current bug that publish doesn't work without an output closure
output {
    "." {
        mode "copy"
    }
    "annotation" {
        mode "copy"
    }
    "variation" {
        mode "copy"
    }
}

