#!/usr/bin/env nextflow

/* 
 * NextFlow test pipe
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 * 
 */

/*
 * Input parameters: read pairs
 * Params are stored in the params.config file
 */

version                 = "1.0"
// this prevents a warning of undefined parameter
params.help             = false

// this prints the input parameters
log.info """
BIOCORE@CRG - N F TESTPIPE  ~  version ${version}
=============================================
reads                           : ${params.reads}
reference                       : ${params.reference}
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
fastqcOutputFolder    = "ouptut_fastqc"
alnOutputFolder       = "ouptut_aln"
multiqcOutputFolder   = "ouptut_multiQC"

 
Channel
    .fromFilePairs( params.reads )  											 // read the files indicated by the wildcard                            
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" } // if empty, complains
    .set {reads} 														 // make the channel "reads_for_fastqc"

reads.map{
		[it[1]]
	}.collect()
     .view()

reference = file(params.reference)
overhang = params.overhang
annotation = file(params.annotation)

include { FASTQC } from "${baseDir}/modules/fastqc" addParams(OUTPUT: fastqcOutputFolder, LABEL:"fourcpus")
include { STAR } from "${baseDir}/modules/star" addParams(OUTPUT: alnOutputFolder, LABEL:"fourcpus")
include { MATRIX } from "${baseDir}/modules/star" addParams(OUTPUT: alnOutputFolder, LABEL:"fourcpus")

workflow {
	fastqc_out = FASTQC(reads)
	map_res = STAR(overhang, reference, annotation, reads)
        map_res.quants.view()
	map_res.quants.collect().view()
        matrix = MATRIX(map_res.quants)
}



workflow.onComplete { 
	println ( workflow.success ? "\nDone! Open the following report in your browser --> ${multiqcOutputFolder}/multiqc_report.html\n" : "Oops .. something went wrong" )
}

