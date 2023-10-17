/*
*  Trimmomatic module
*/

params.CONTAINER = "quay.io/biocontainers/trimmomatic:0.32--hdfd78af_4"
params.OUTPUT = "trimmomatic_output"
params.LABEL = ""
params.EXTRAPARS = ""

process trim_PE {
    publishDir(params.OUTPUT, mode: 'copy')
    tag {"${pair_id}"}
    container params.CONTAINER
    label params.LABEL

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}_*_paired-trimmed.fastq.gz"), emit: trimmed_reads //there is only one ouput, emit not necessary

    """
    trimmomatic PE -threads ${task.cpus} \
    -trimlog trim.logs \
    ${reads} \
    ${pair_id}_1_paired-trimmed.fastq.gz ${pair_id}_1_unpaired-trimmed.fastq.gz \
    ${pair_id}_2_paired-trimmed.fastq.gz ${pair_id}_2_unpaired-trimmed.fastq.gz \
    ${params.EXTRAPARS}
    """
}

workflow TRIM_PE{
    take:
    reads

    main:
    trim_pe = trim_PE(reads)
    
    emit:
    trimpe = trim_pe.trimmed_reads
}

