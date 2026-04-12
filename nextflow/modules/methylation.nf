// modules/methylation.nf — modkit pileup (5mCG/5hmCG)

process MODKIT_PILEUP {
    tag "${sample}"
    publishDir "${params.outdir}/methylation", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path ref

    output:
    tuple val(sample), path("${sample}.bedmethyl.gz"),     emit: bedmethyl
    tuple val(sample), path("${sample}.bedmethyl.gz.tbi"), emit: tbi
    path "${sample}.modkit_summary.txt",                   emit: summary
    path "${sample}.modkit.log",                           emit: log

    script:
    """
    # Summary of modification calls
    modkit summary \
        ${bam} \
        --threads ${task.cpus} \
        --log ${sample}.modkit.log \
        > ${sample}.modkit_summary.txt 2>> ${sample}.modkit.log || true

    # Per-CpG methylation pileup (5mCG + 5hmCG)
    modkit pileup \
        ${bam} \
        ${sample}.bedmethyl \
        --ref ${ref} \
        --modified-bases 5mC 5hmC \
        --cpg \
        --threads ${task.cpus} \
        --log ${sample}.modkit.log \
        --filter-threshold 0.66

    # Compress and index
    bgzip ${sample}.bedmethyl
    tabix -p bed ${sample}.bedmethyl.gz

    echo "modkit pileup complete for ${sample}"
    zcat ${sample}.bedmethyl.gz | head -5
    """
}
