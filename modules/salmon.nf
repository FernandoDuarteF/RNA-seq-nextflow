/*
* Salmon module
*/

params.CONTAINER = "quay.io/biocontainers/salmon:1.10.2--hecfa306_0ii"
params.OUTPUT = "salmon_output"
params.LABEL = ""

process lib_type {
    publishDir(params.OUTPUT, mode: 'copy')
    tag "${ref.getSimpleName()}"
    container params.CONTAINER
    label params.LABEL

    input:
    val overhang
    path ref
    path anno

    output:
    path "${ref.getSimpleName()}" //there is only one ouput, emit not necessary

    """
    mkdir ${ref.getSimpleName()}
    STAR --runThreadN ${task.cpus} \
    --runMode genomeGenerate \
    --genomeDir ${ref.getSimpleName()} \
    --genomeFastaFiles ${ref} \
    --sjdbOverhang ${overhang} \
    --sjdbGTFfile ${anno} \
    --genomeSAindexNbases 11
    """
}

process lib_type {
    publishDir(params.OUTPUT, mode: 'copy')
    tag "lib_type"
    container params.CONTAINER
    label params.LABEL

    input:
    val overhang
    path ref
    path anno

    output:
    path "${ref.getSimpleName()}" //there is only one ouput, emit not necessary

    """
    mkdir ${ref.getSimpleName()}
    STAR --runThreadN ${task.cpus} \
    --runMode genomeGenerate \
    --genomeDir ${ref.getSimpleName()} \
    --genomeFastaFiles ${ref} \
    --sjdbOverhang ${overhang} \
    --sjdbGTFfile ${anno} \
    --genomeSAindexNbases 11
    """
}
