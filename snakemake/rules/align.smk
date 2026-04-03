# rules/align.smk — Minimap2 + samtools sort/index

rule minimap2_align:
    input:
        fastq = lambda wc: config["samples"][wc.sample]["fastq"],
        ref   = config["ref"],
    output:
        sam = temp("{outdir}/alignment/{sample}.sam"),
    params:
        preset = config["minimap2"]["preset"],
        extra  = config["minimap2"]["extra"],
        rg     = "@RG\\tID:{sample}\\tSM:{sample}\\tPL:ONT",
    threads: config["threads"]
    conda:  "envs/align.yaml"
    log:    "{outdir}/logs/minimap2_{sample}.log"
    shell:
        """
        minimap2 \
            -ax {params.preset} \
            -t {threads} \
            {params.extra} \
            -R "{params.rg}" \
            {input.ref} \
            {input.fastq} \
            > {output.sam} \
            2> {log}
        """

rule samtools_sort_index:
    input:
        # Use existing BAM if provided, else use minimap2 output
        bam = lambda wc: (
            config["samples"][wc.sample]["bam"]
            if config["samples"][wc.sample].get("bam")
            else f"{OUTDIR}/alignment/{wc.sample}.sam"
        )
    output:
        bam   = "{outdir}/alignment/{sample}.sorted.bam",
        bai   = "{outdir}/alignment/{sample}.sorted.bam.bai",
        stats = "{outdir}/alignment/{sample}.flagstat.txt",
        idx   = "{outdir}/alignment/{sample}.idxstats.txt",
    threads: config["threads"]
    conda:  "envs/align.yaml"
    log:    "{outdir}/logs/samtools_{sample}.log"
    shell:
        """
        # Sort (accepts SAM or BAM)
        samtools sort \
            -@ {threads} \
            -o {output.bam} \
            {input.bam} \
            2>> {log}

        # Index
        samtools index -@ {threads} {output.bam} 2>> {log}

        # Statistics
        samtools flagstat {output.bam} > {output.stats} 2>> {log}
        samtools idxstats {output.bam} > {output.idx}   2>> {log}

        echo "Flagstat for {wildcards.sample}:" >> {log}
        cat {output.stats} >> {log}
        """
