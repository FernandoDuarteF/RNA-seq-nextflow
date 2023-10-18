/*
*  STAR module
*/

params.CONTAINER = "quay.io/biocontainers/star:2.7.7a--0"
params.OUTPUT = "star_output"
params.LABEL = ""
params.EXTRAPARS = ""

process index {
    publishDir(params.OUTPUT, mode: 'copy')
    tag {"${ref.getSimpleName()}"}
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
    --sjdbGTFfile ${anno}
    """
}

process aln {
    publishDir(params.OUTPUT, mode: 'copy', pattern: '*ReadsPerGene.out.tab')
    tag {"${pair_id}"}
    container params.CONTAINER
    label params.LABEL

    input:
    path index
    tuple val(pair_id), path (reads)

    output:
    tuple val(pair_id), path("${pair_id}*.bam"), emit: bams 
    tuple val(pair_id), path("${pair_id}ReadsPerGene.out.tab"), emit: quants
    tuple val(pair_id), path("${pair_id}Log.final.out"), emit: logs 
    tuple val(pair_id), path("${pair_id}SJ*"),  emit: junctions
    // use \ to scape $ and use them in bash
    script:
    """
    commandgz="" && [ \$(echo "${reads}" | cut -d" " -f1 | grep ".gz") ] && commandgz="--readFilesCommand zcat "
    STAR --genomeDir ${index} \
    --readFilesIn ${reads} \
    \$commandgz \
    --runThreadN ${task.cpus} \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --outFileNamePrefix ${pair_id} \
    ${params.EXTRAPARS}
    """

}

process matrix {

    publishDir(params.OUTPUT, mode: 'copy')
    tag "matrix"
    label params.LABEL

    input:
    val col
    path counts

    output:
    path "raw_counts_matrix.txt"

    script:
    """
    paste ${counts} | grep -v "_" | awk '{printf "%s\t", \$1}{for (i=${col};i<=NF;i+=${col}) printf "%s\t", \$i; printf "\\n" }' > tmp
    sed -e "1igene_name\t\$(ls ${counts} | tr '\n' '\t' | sed 's/ReadsPerGene.out.tab//g')" tmp | cut -f1-9 > raw_counts_matrix.txt
    rm tmp
    """

}


workflow STAR {
    take:
    overhang
    reference
    annotation
    reads

    main:
    star_index = index(overhang, reference, annotation)
    star_aln = aln(star_index, reads)

    emit:
    quants = star_aln.quants
    logs = star_aln.logs
}

workflow MATRIX {
    take:
    counts

    main:
    counts.map{
                [it[1]] //Remove id the from emited tuple and collect output counts
        }.collect().set{countsm}
    out = matrix(countsm)
    
    emit:
    out
}
