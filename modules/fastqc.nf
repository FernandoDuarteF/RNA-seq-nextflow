/*
*  fastqc module
*/

params.CONTAINER = "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
params.OUTPUT = "fastq_output"
params.LABEL = ""

process fastqc {
    publishDir(params.OUTPUT, mode: 'copy')
    tag { "${reads}" }
    container params.CONTAINER
    label (params.LABEL)

    input:
    path(reads)

    output:
    path("*_fastqc*")

    script:
    """
        fastqc -t ${task.cpus} ${reads}
    """
}

workflow FASTQC {
    take:
    reads

    main:
    reads.map{
		[it[1]] //fromFailePairs output a tuple with id and both pair reads. Discard the id
	}.flatten().set{readsqc} //flatten tuple, equivalent to extend all fastq files (check .view() if not sure)
    fastqc(readsqc)
}