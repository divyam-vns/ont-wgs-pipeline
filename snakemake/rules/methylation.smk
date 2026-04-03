# rules/methylation.smk — modkit pileup (5mCG/5hmCG)

rule modkit_summary:
    input:
        bam = "{outdir}/alignment/{sample}.sorted.bam",
        bai = "{outdir}/alignment/{sample}.sorted.bam.bai",
    output:
        summary = "{outdir}/methylation/{sample}.modkit_summary.txt",
    threads: config["threads"]
    conda:  "envs/methylation.yaml"
    log:    "{outdir}/logs/modkit_summary_{sample}.log"
    shell:
        """
        modkit summary \
            {input.bam} \
            --threads {threads} \
            --log-filepath {log} \
            > {output.summary}
        """

rule modkit_pileup:
    input:
        bam     = "{outdir}/alignment/{sample}.sorted.bam",
        bai     = "{outdir}/alignment/{sample}.sorted.bam.bai",
        ref     = config["ref"],
        summary = "{outdir}/methylation/{sample}.modkit_summary.txt",
    output:
        bedmethyl     = "{outdir}/methylation/{sample}.bedmethyl",
        bedmethyl_gz  = "{outdir}/methylation/{sample}.bedmethyl.gz",
        bedmethyl_tbi = "{outdir}/methylation/{sample}.bedmethyl.gz.tbi",
    params:
        filter_threshold = config["modkit"]["filter_threshold"],
        preset           = config["modkit"]["preset"],
        extra            = config["modkit"]["extra"],
    threads: config["threads"]
    conda:  "envs/methylation.yaml"
    log:    "{outdir}/logs/modkit_pileup_{sample}.log"
    shell:
        """
        modkit pileup \
            {input.bam} \
            {output.bedmethyl} \
            --ref {input.ref} \
            --threads {threads} \
            --log-filepath {log} \
            --preset {params.preset} \
            --filter-threshold {params.filter_threshold} \
            --mod-thresholds m:{params.filter_threshold} \
            {params.extra}

        # Compress + tabix index
        bgzip {output.bedmethyl}
        tabix -p bed {output.bedmethyl_gz}

        echo "Top 5 CpG methylation calls:"
        zcat {output.bedmethyl_gz} | head -5
        """
