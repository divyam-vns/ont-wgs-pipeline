# rules/variant_calling.smk — Clair3 (SNP/indel) + Sniffles2 (SV)

# ---- Clair3 SNP/indel ----
rule clair3_snp:
    input:
        bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai = "{outdir}/alignment/{sample}.sorted.bam.bai",
        ref = config["ref"],
    output:
        vcf = "{outdir}/variants/snp/{sample}/merge_output.vcf.gz",
        tbi = "{outdir}/variants/snp/{sample}/merge_output.vcf.gz.tbi",
    params:
        outdir      = "{outdir}/variants/snp/{sample}",
        platform    = config["clair3"]["platform"],
        model_path  = config["clair3"]["model_path"],
        extra       = config["clair3"]["extra"],
        bed_arg     = f"--bed_fn={config['bed']}" if config.get("bed") else "",
    threads: config["threads"]
    conda:  "envs/variants.yaml"
    log:    "{outdir}/logs/clair3_{sample}.log"
    shell:
        """
        # Ensure reference is indexed
        if [ ! -f {input.ref}.fai ]; then
            samtools faidx {input.ref}
        fi

        run_clair3.sh \
            --bam_fn={input.bam} \
            --ref_fn={input.ref} \
            --threads={threads} \
            --platform={params.platform} \
            --model_path={params.model_path} \
            --output={params.outdir} \
            --sample_name={wildcards.sample} \
            {params.bed_arg} \
            {params.extra} \
            > {log} 2>&1
        """

# ---- Sniffles2 SV ----
rule sniffles2_sv:
    input:
        bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai = "{outdir}/alignment/{sample}.sorted.bam.bai",
        ref = config["ref"],
    output:
        vcf     = "{outdir}/variants/sv/{sample}.sniffles.vcf",
        vcf_gz  = "{outdir}/variants/sv/{sample}.sniffles.vcf.gz",
        snf     = "{outdir}/variants/sv/{sample}.sniffles.snf",
    params:
        min_sv_len = config["sniffles2"]["min_sv_len"],
        min_mapq   = config["sniffles2"]["min_mapq"],
        extra      = config["sniffles2"]["extra"],
    threads: config["threads"]
    conda:  "envs/variants.yaml"
    log:    "{outdir}/logs/sniffles2_{sample}.log"
    shell:
        """
        sniffles \
            --input {input.bam} \
            --vcf {output.vcf} \
            --snf {output.snf} \
            --reference {input.ref} \
            --threads {threads} \
            --sample-id {wildcards.sample} \
            --minsvlen {params.min_sv_len} \
            --mapq {params.min_mapq} \
            {params.extra} \
            > {log} 2>&1

        # Compress + index
        bgzip -c {output.vcf} > {output.vcf_gz}
        tabix -p vcf {output.vcf_gz}
        """
