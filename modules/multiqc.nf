/*
*  multiqc module
*/

params.CONTAINER = "quay.io/biocontainers/multiqc:1.17--pyhdfd78af_0"
params.OUTPUT = "multiqc_output"
params.LABEL = ""

process multiqc {
    publishDir(params.OUTPUT, mode: 'copy')
    tag "multiqc"
    container params.CONTAINER
    label params.LABEL

    input:
    path (inputfolders)

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc -f ${inputfolders}
    """
}
