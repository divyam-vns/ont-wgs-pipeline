// modules/align.nf — Minimap2 alignment + samtools sort/index

process MINIMAP2_ALIGN {
    tag "${sample}"
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: '*.sam'

    input:
    tuple val(sample), path(reads)
    path ref

    output:
    tuple val(sample), path("${sample}.sam"), emit: sam

    script:
    """
    minimap2 \
        -ax map-ont \
        -t ${task.cpus} \
        --MD \
        --secondary=no \
        -R "@RG\\tID:${sample}\\tSM:${sample}\\tPL:ONT" \
        ${ref} \
        ${reads} \
        > ${sample}.sam

    echo "Alignment complete: ${sample}.sam"
    """
}

process SAMTOOLS_SORT_INDEX {
    tag "${sample}"
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(sample), path(sam)

    output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: bam

    script:
    """
    samtools sort \
        -@ ${task.cpus} \
        -o ${sample}.sorted.bam \
        ${sam}

    samtools index \
        -@ ${task.cpus} \
        ${sample}.sorted.bam

    samtools flagstat ${sample}.sorted.bam > ${sample}.flagstat.txt
    samtools idxstats ${sample}.sorted.bam > ${sample}.idxstats.txt

    echo "Flagstat:"
    cat ${sample}.flagstat.txt
    """
}
