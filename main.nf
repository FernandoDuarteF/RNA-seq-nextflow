#!/usr/bin/env nextflow

/* 
 * Author: Fernando Duarte
 * Email : <d.fernando.fg@gmail.com>
 */

version                 = "1.0"
// this prevents a warning of undefined parameter
params.help             = false

// this prints the input parameters
log.info """
RNA-seq pipeline  ~  version ${version}
=============================================
reads                           : ${params.reads}
reference                       : ${params.reference}
annotation                      : ${params.annotation}
overhang                        : ${params.overhang}
"""

// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
    log.info 'This is the Biocore\'s NF test pipeline'
    log.info 'Enjoy!'
    log.info '\n'
    exit 1
}

/*
 * Defining the output folders.
 */
fastqcOutputFolder    = "results/fastqc"
alnOutputFolder       = "results/aln"
multiqcOutputFolder   = "results/multiQC"
trimOutputFolder      = "results/trimmomatic"
 
Channel
    .fromFilePairs( params.reads )  											 // read the files indicated by the wildcard                            
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" } // if empty, complains
    .set {reads} 														 // make the channel "reads"

// to check
reads.map{
		[it[1]]
	}.collect()
     .view()

reference  = file(params.reference)
overhang   = params.overhang
annotation = file(params.annotation)
adaptors   = file(params.adaptors)

include { FASTQC } from "${baseDir}/modules/fastqc" addParams(OUTPUT: fastqcOutputFolder, LABEL:"qc")
include { FASTQC as FASTQC_TRIM } from "${baseDir}/modules/fastqc" addParams(OUTPUT: fastqcOutputFolder, LABEL:"qc")
include { TRIM_PE } from "${baseDir}/modules/trimmomatic.nf" addParams(OUTPUT: trimOutputFolder, LABEL:"idx", EXTRAPARS:"ILLUMINACLIP:${adaptors}:2:30:10")
include { STAR } from "${baseDir}/modules/star" addParams(OUTPUT: alnOutputFolder, LABEL:"idx")
include { MATRIX } from "${baseDir}/modules/star" addParams(OUTPUT: alnOutputFolder)
include { multiqc } from "${baseDir}/modules/multiqc" addParams(OUTPUT: multiqcOutputFolder)

workflow {
	fastqc_out = FASTQC(reads)
  trimmomatic = TRIM_PE(reads)
  reads_trim = trimmomatic.trimpe
  reads_trim_qc = FASTQC_TRIM(reads_trim)
  //reads_trim.view()
	map_res = STAR(overhang, reference, annotation, reads_trim)
  to_multiqc = reads_trim_qc.qc_out.mix(map_res.quants.map{[it[1]]}).mix(map_res.logs.map{[it[1]]}).collect()
  to_multiqc.view()
  multiqc(to_multiqc)
  //map_res.quants.view()
	//map_res.quants.collect().view()
       // matrix = MATRIX(map_res.quants)
}



workflow.onComplete { 
	println ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
}

