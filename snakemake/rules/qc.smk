# rules/qc.smk — NanoPlot + MultiQC

rule nanoplot:
    input:
        bam = lambda wc: config["samples"][wc.sample]["bam"]
    output:
        stats = "{outdir}/qc/{sample}/NanoStats.txt",
        html  = "{outdir}/qc/{sample}/NanoPlot-report.html",
    params:
        outdir  = "{outdir}/qc/{sample}",
        extra   = config["nanoplot"]["extra"],
    threads: config["threads"]
    conda:  "envs/align.yaml"
    log:    "{outdir}/logs/nanoplot_{sample}.log"
    shell:
        """
        NanoPlot \
            --bam {input.bam} \
            --sample_name {wildcards.sample} \
            --threads {threads} \
            --outdir {params.outdir} \
            {params.extra} \
            > {log} 2>&1
        """

rule multiqc:
    input:
        expand("{outdir}/qc/{sample}/NanoStats.txt", sample=SAMPLES, outdir=OUTDIR)
    output:
        report = "{outdir}/qc/multiqc_report.html",
        data   = directory("{outdir}/qc/multiqc_data"),
    params:
        qc_dir = "{outdir}/qc",
    conda:  "envs/align.yaml"
    log:    "{outdir}/logs/multiqc.log"
    shell:
        """
        multiqc {params.qc_dir} \
            --title "ONT WGS QC Report" \
            --outdir {params.qc_dir} \
            --filename multiqc_report.html \
            --force \
            > {log} 2>&1
        """
