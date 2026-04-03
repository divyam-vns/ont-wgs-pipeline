// modules/qc.nf — NanoPlot + MultiQC

process NANOPLOT {
    tag "${sample}"
    publishDir "${params.outdir}/qc/nanoplot/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    path "*.html",              emit: html
    path "NanoStats.txt",       emit: stats
    path "*.png",   optional: true

    script:
    def input_flag = reads.name.endsWith('.bam') ? "--bam ${reads}" : "--fastq ${reads}"
    """
    NanoPlot \
        ${input_flag} \
        --sample_name ${sample} \
        --threads ${task.cpus} \
        --outdir . \
        --plots kde dot \
        --N50 \
        --loglength \
        --no_static
    """
}

process MULTIQC {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path stats_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/",       emit: data

    script:
    """
    multiqc . \
        --title "ONT WGS Pipeline QC" \
        --filename multiqc_report.html \
        --force
    """
}
