// modules/methylation.nf — modkit pileup for 5mCG/5hmCG detection

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
    # Check BAM has MM/ML modification tags
    MOD_TAGS=\$(samtools view -H ${bam} | grep -c "MM:Z" || true)
    if [ "\$MOD_TAGS" -eq 0 ]; then
        echo "WARNING: BAM may not contain MM/ML modification tags. Proceeding anyway..."
    fi

    # Summary of modification calls
    modkit summary \
        ${bam} \
        --threads ${task.cpus} \
        --log-filepath ${sample}.modkit.log \
        > ${sample}.modkit_summary.txt

    # Per-CpG methylation pileup
    modkit pileup \
        ${bam} \
        ${sample}.bedmethyl \
        --ref ${ref} \
        --threads ${task.cpus} \
        --log-filepath ${sample}.modkit.log \
        --preset traditional \
        --filter-threshold 0.66 \
        --mod-thresholds m:0.66 \
        --cpg \
        --ignore h

    # Compress and index
    bgzip ${sample}.bedmethyl
    tabix -p bed ${sample}.bedmethyl.gz

    echo "modkit pileup complete for ${sample}"
    echo "Top 10 CpG sites:"
    zcat ${sample}.bedmethyl.gz | head -10
    """
}
