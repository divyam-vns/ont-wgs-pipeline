# Pipeline Diagram

## Full Analysis DAG

```
ONT Demo BAM (R10.4.1 HAC, MM/ML tags)
              │
    ┌─────────┴──────────┐
    │                    │
  [QC]               [Already aligned?]
  NanoPlot            Yes → pass through
  MultiQC             No  → Minimap2 (-ax map-ont)
    │                         │
    │                   samtools sort/index
    │                         │
    └──────────┬──────────────┘
               │
         sorted BAM + BAI
               │
    ┌──────────┼──────────────┐
    │          │              │
 [SNP/indel] [SV]       [Methylation]
  Clair3    Sniffles2    modkit pileup
    │          │              │
  VCF.gz   VCF.gz        bedMethyl.gz
  TBI        TBI            TBI
    │          │              │
    └──────────┴──────────────┘
                    │
             summarize_results.py
                    │
            summary TSV + report
```

## Module Details

### QC Module
```
Input: BAM or FASTQ
  ↓
NanoPlot --bam / --fastq
  → NanoStats.txt (N50, median length, median quality, total bases)
  → NanoPlot-report.html
  ↓
MultiQC (aggregates all NanoStats)
  → multiqc_report.html
```

### Alignment Module (FASTQ input only)
```
Input: FASTQ + Reference FASTA
  ↓
minimap2 -ax map-ont --MD --secondary=no
  → {sample}.sam
  ↓
samtools sort -@ threads
  → {sample}.sorted.bam
  ↓
samtools index
  → {sample}.sorted.bam.bai
  ↓
samtools flagstat / idxstats
  → {sample}.flagstat.txt
```

### SNP/Indel Module
```
Input: sorted BAM + Reference FASTA + (optional BED)
  ↓
Clair3 run_clair3.sh
  --platform ont
  --model_path dna_r10.4.1_e8.2_400bps_hac@v4.2.0
  → pileup.vcf.gz (pileup model)
  → full_alignment.vcf.gz
  ↓
merge_output.vcf.gz (final merged calls)
  → merge_output.vcf.gz.tbi
```

### Structural Variant Module
```
Input: sorted BAM + Reference FASTA
  ↓
Sniffles2
  --minsvlen 50 --mapq 20 --phase --output-rnames
  → {sample}.sniffles.vcf (DEL, INS, DUP, INV, BND)
  → {sample}.sniffles.snf (for multi-sample merging)
  ↓
bgzip + tabix
  → {sample}.sniffles.vcf.gz + .tbi
```

### Methylation Module
```
Input: sorted BAM with MM/ML tags + Reference FASTA
  ↓
modkit summary
  → {sample}.modkit_summary.txt (base modification overview)
  ↓
modkit pileup
  --preset traditional
  --cpg --ignore h
  --filter-threshold 0.66
  → {sample}.bedmethyl (per-CpG: chrom, start, end, mod_code,
                         n_valid_cov, fraction_modified)
  ↓
bgzip + tabix
  → {sample}.bedmethyl.gz + .tbi
```

## Key Tool Versions

| Tool      | Version | Container                                          |
|-----------|---------|----------------------------------------------------|
| minimap2  | ≥2.26   | quay.io/biocontainers/minimap2:2.26--he4a0461_2    |
| samtools  | 1.23    | (same as minimap2 container)                       |
| Clair3    | ≥1.0.4  | docker.io/hkubal/clair3:latest                     |
| Sniffles2 | ≥2.3    | quay.io/biocontainers/sniffles:2.3.3--pyhdfd78af_0 |
| modkit    | ≥0.3    | docker.io/ontresearch/modkit:latest                |
| NanoPlot  | ≥1.42   | quay.io/biocontainers/nanoplot:1.42.0--pyhdfd78af_0|

## Demo Data Chemistry

| Property         | Value                                      |
|------------------|--------------------------------------------|
| Source           | ONT EPI2ME Labs wf-human-variation demo    |
| Chemistry        | R10.4.1 e8.2 400bps                        |
| Basecaller model | dna_r10.4.1_e8.2_400bps_hac@v4.2.0        |
| Mod model        | _5mCG_5hmCG@v2                             |
| Reference        | GRCh38 (subset)                            |
| MM/ML tags       | Yes (5mCG + 5hmCG)                         |
