// modules/variant_calling.nf — Clair3 (SNP/indel) + Sniffles2 (SV)

process CLAIR3_SNP {
    tag "${sample}"
    publishDir "${params.outdir}/variants/snp", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path ref
    path bed   // optional; use [] if not provided

    output:
    tuple val(sample), path("merge_output.vcf.gz"),     emit: vcf
    tuple val(sample), path("merge_output.vcf.gz.tbi"), emit: tbi
    path "run_clair3.log",                              emit: log

    script:
    def bed_arg = bed ? "--bed_fn=${bed}" : ""
    """
    # Index reference if needed
    if [ ! -f ${ref}.fai ]; then
        samtools faidx ${ref}
    fi

    # Find built-in model path inside container
    MODEL=\$(find /opt/models -maxdepth 1 -type d | grep -v "^/opt/models\$" | head -1)
    if [ -z "\$MODEL" ]; then
        MODEL=\$(find /usr -name "*.cfg" -path "*/clair3_models/*" 2>/dev/null | head -1 | xargs dirname || echo "")
    fi
    if [ -z "\$MODEL" ]; then
        MODEL=\$(ls /opt/conda/bin/../share/clair3*/models/ 2>/dev/null | head -1 || echo "")
        MODEL="/opt/conda/share/clair3/models/\$MODEL"
    fi
    echo "Using model: \$MODEL"

    run_clair3.sh \
        --bam_fn=${bam} \
        --ref_fn=${ref} \
        --threads=${task.cpus} \
        --platform=ont \
        --model_path=\$MODEL \
        --output=. \
        --sample_name=${sample} \
        ${bed_arg} \
        --include_all_ctgs \
        --no_phasing_for_fa \
        2>&1 | tee run_clair3.log

    echo "Clair3 SNP/indel calling complete for ${sample}"
    """
}

process SNIFFLES2_SV {
    tag "${sample}"
    publishDir "${params.outdir}/variants/sv", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path ref

    output:
    tuple val(sample), path("${sample}.sniffles.vcf"),    emit: vcf
    tuple val(sample), path("${sample}.sniffles.snf"),    emit: snf
    path "${sample}.sniffles.log",                         emit: log

    script:
    """
    sniffles \
        --input ${bam} \
        --vcf ${sample}.sniffles.vcf \
        --snf ${sample}.sniffles.snf \
        --reference ${ref} \
        --threads ${task.cpus} \
        --sample-id ${sample} \
        --minsvlen 50 \
        --mapq 20 \
        --output-rnames \
        --phase \
        2>&1 | tee ${sample}.sniffles.log

    # Compress and index
    bgzip -c ${sample}.sniffles.vcf > ${sample}.sniffles.vcf.gz
    tabix -p vcf ${sample}.sniffles.vcf.gz

    echo "Sniffles2 SV calling complete for ${sample}"
    """
}
