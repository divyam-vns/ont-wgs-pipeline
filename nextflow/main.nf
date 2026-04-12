#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================
// ONT WGS Pipeline — Nextflow DSL2
// Modules: QC → Alignment → Variant Calling → Methylation
// ============================================================

// --------------- Parameters ---------------
params.input         = null          // BAM or FASTQ input
params.ref           = null          // Reference FASTA (GRCh38 subset for demo)
params.bed           = null          // Optional target BED
params.sample        = "SAMPLE"
params.outdir        = "results/nextflow"
params.basecaller_cfg = "dna_r10.4.1_e8.2_400bps_hac@v4.2.0"
params.mod_cfg       = "dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2"
params.threads       = 8
params.min_coverage  = 5

// Sub-workflow flags
params.run_qc        = true
params.run_snp       = true
params.run_sv        = true
params.run_methylation = true

// --------------- Imports ---------------
include { NANOPLOT; MULTIQC }          from './modules/qc.nf'
include { MINIMAP2_ALIGN; SAMTOOLS_SORT_INDEX } from './modules/align.nf'
include { CLAIR3_SNP; SNIFFLES2_SV }   from './modules/variant_calling.nf'
include { MODKIT_PILEUP }              from './modules/methylation.nf'

// --------------- Input validation ---------------
def check_params() {
    if (!params.input) error "ERROR: --input is required (BAM or FASTQ)"
    if (!params.ref)   error "ERROR: --ref is required (reference FASTA)"
}

// --------------- Workflow ---------------
workflow {
    check_params()

    // Determine if input is BAM or FASTQ
    input_file  = file(params.input)
    ref_file    = file(params.ref)
    bed_file    = params.bed ? file(params.bed) : []

    // ---- QC ----
    if (params.run_qc) {
        NANOPLOT(
            Channel.of( [params.sample, input_file] )
        )
        MULTIQC(
            NANOPLOT.out.stats.collect()
        )
    }

    // ---- Alignment (only if FASTQ input) ----
    def bam_ch
    if (input_file.name.endsWith('.bam') || input_file.name.endsWith('.cram')) {
        // Already aligned — find BAI and pass [sample, bam, bai] tuple
        def bai_file = file(params.input + '.bai')
        bam_ch = Channel.of( [params.sample, input_file, bai_file] )
    } else {
        MINIMAP2_ALIGN(
            Channel.of( [params.sample, input_file] ),
            ref_file
        )
        SAMTOOLS_SORT_INDEX( MINIMAP2_ALIGN.out.sam )
        bam_ch = SAMTOOLS_SORT_INDEX.out.bam
    }

    // ---- Variant Calling ----
    if (params.run_snp) {
        CLAIR3_SNP(
            bam_ch,
            ref_file,
            bed_file
        )
    }

    if (params.run_sv) {
        SNIFFLES2_SV(
            bam_ch,
            ref_file
        )
    }

    // ---- Methylation ----
    if (params.run_methylation) {
        MODKIT_PILEUP(
            bam_ch,
            ref_file
        )
    }

    // ---- Summary ----
    workflow.onComplete {
        log.info """
        ============================================
        ONT WGS Pipeline — COMPLETE
        ============================================
        Sample   : ${params.sample}
        Input    : ${params.input}
        Reference: ${params.ref}
        Results  : ${params.outdir}
        Duration : ${workflow.duration}
        Status   : ${workflow.success ? 'SUCCESS ✓' : 'FAILED ✗'}
        ============================================
        """.stripIndent()
    }
}
