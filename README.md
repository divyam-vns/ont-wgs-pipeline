# ONT WGS Pipeline 🧬

**Human Whole Genome Sequencing analysis using Oxford Nanopore Technology**  
Dual-framework implementation: **Nextflow (DSL2)** + **Snakemake**

---

## Overview

This pipeline processes ONT long-read WGS data through four analysis modules:

| Module | Tools | Output |
|--------|-------|--------|
| **QC** | NanoPlot, FastQC, MultiQC | HTML reports, quality metrics |
| **Alignment** | Minimap2 (map-ont), samtools | Sorted, indexed BAM |
| **Variant Calling** | Clair3 (SNP/indel), Sniffles2 (SV) | VCF/BCF files |
| **Methylation** | modkit | bedMethyl, CpG methylation profiles |

### Demo Dataset
Uses the **ONT open data wf-human-variation demo** (~small chr subset, R10.4.1 chemistry):
```bash
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-human-variation/demo_data.tar.gz
```
Sample: HG002-like demo BAM with MM/ML methylation tags, GRCh38 subset.

---

## Repository Structure

```
ont-wgs-pipeline/
├── nextflow/
│   ├── main.nf                  # Nextflow DSL2 entry point
│   ├── nextflow.config          # Resource profiles (standard, docker, hpc)
│   ├── modules/
│   │   ├── qc.nf                # NanoPlot + MultiQC
│   │   ├── align.nf             # Minimap2 alignment
│   │   ├── variant_calling.nf   # Clair3 + Sniffles2
│   │   └── methylation.nf       # modkit pileup
│   └── conf/
│       ├── docker.config
│       └── resources.config
├── snakemake/
│   ├── Snakefile                # Main Snakemake workflow
│   ├── config.yaml              # Sample sheet + parameters
│   ├── rules/
│   │   ├── qc.smk
│   │   ├── align.smk
│   │   ├── variant_calling.smk
│   │   └── methylation.smk
│   └── envs/
│       ├── align.yaml
│       ├── variants.yaml
│       └── methylation.yaml
├── scripts/
│   ├── download_demo.sh         # Download ONT demo data + reference
│   ├── setup_env.sh             # Conda/Podman environment setup
│   └── summarize_results.py     # Aggregate QC + variant stats
├── docs/
│   └── pipeline_diagram.md
└── data/demo/                   # Demo data lands here (gitignored)
```

---

## Quick Start

### 1. Setup Environment
```bash
git clone https://github.com/divyam-vns/ont-wgs-pipeline.git
cd ont-wgs-pipeline
bash scripts/setup_env.sh
```

### 2. Download Demo Data
```bash
bash scripts/download_demo.sh
# Downloads ~300 MB demo BAM + reference subset into data/demo/
```

### 3. Run — Nextflow
```bash
# With Docker (recommended)
nextflow run nextflow/main.nf \
  --input   data/demo/demo.bam \
  --ref     data/demo/demo.fasta \
  --bed     data/demo/demo.bed \
  --sample  DEMO \
  --outdir  results/nextflow \
  -profile  docker

# With Podman (Mac/Linux, no Docker daemon required)
nextflow run nextflow/main.nf \
  --input   data/demo/demo.bam \
  --ref     data/demo/demo.fasta \
  --bed     data/demo/demo.bed \
  --sample  DEMO \
  --outdir  results/nextflow \
  -profile  podman
```

### 4. Run — Snakemake
```bash
# Dry run first
snakemake -s snakemake/Snakefile \
  --configfile snakemake/config.yaml \
  --cores 8 \
  --use-conda \
  -n

# Execute
snakemake -s snakemake/Snakefile \
  --configfile snakemake/config.yaml \
  --cores 8 \
  --use-conda
```

---

## Requirements

| Tool | Version | Install |
|------|---------|---------|
| Nextflow | ≥23.10 | `curl -s get.nextflow.io | bash` |
| Snakemake | ≥7.32 | `conda install -c bioconda snakemake` |
| Docker / Podman | Latest | See `scripts/setup_env.sh` |
| minimap2 | ≥2.26 | bioconda |
| samtools | 1.23 | bioconda |
| Clair3 | ≥1.0 | Docker: `hkubal/clair3:latest` |
| Sniffles2 | ≥2.3 | bioconda |
| modkit | ≥0.3 | ONT / bioconda |
| NanoPlot | ≥1.42 | bioconda |

---

## Pipeline Details

### QC
- **NanoPlot** on input FASTQ/BAM: read length distribution, quality histogram, N50
- **FastQC** for secondary QC
- **MultiQC** aggregates all reports into a single HTML dashboard

### Alignment
- **Minimap2** with `-ax map-ont` preset for R9/R10 long reads
- MD tags enabled (`--MD`) for downstream variant calling
- BAM sorted + indexed with **samtools**

### Variant Calling
- **Clair3**: germline SNPs and indels; model auto-detected from BAM header  
  (`dna_r10.4.1_e8.2_400bps_hac@v4.2.0` for demo data)
- **Sniffles2**: structural variants (DEL, INS, DUP, INV, BND) ≥50 bp

### Methylation (5mCG/5hmCG)
- **modkit pileup**: per-CpG methylation fraction from MM/ML BAM tags
- Output: bedMethyl format (chrom, start, end, mod_code, n_valid, fraction)
- Demo BAM includes `--remora_cfg dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2` tags

---

## Output Files

```
results/
├── qc/
│   ├── nanoplot/           # NanoPlot HTML + stats
│   └── multiqc_report.html
├── alignment/
│   ├── {sample}.sorted.bam
│   └── {sample}.sorted.bam.bai
├── variants/
│   ├── snp/
│   │   ├── merge_output.vcf.gz      # Clair3 SNP/indel calls
│   │   └── merge_output.vcf.gz.tbi
│   └── sv/
│       ├── {sample}.sniffles.vcf    # Sniffles2 SV calls
│       └── {sample}.sniffles.vcf.gz
└── methylation/
    ├── {sample}.bedmethyl.gz        # Per-CpG methylation
    └── {sample}.bedmethyl.gz.tbi
```

---

## Podman Setup (macOS)

This project is optimized for the **Podman + Miniconda3** stack:
```bash
# Start Podman machine if not running
podman machine start

# Nextflow uses podman automatically with -profile podman
# All containers are pulled from quay.io/biocontainers or Docker Hub
```

---

## References

- [ONT wf-human-variation](https://github.com/epi2me-labs/wf-human-variation)
- [ONT Open Data (AWS)](https://registry.opendata.aws/ont-open-data/)
- [Clair3](https://github.com/HKU-BAL/Clair3)
- [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)
- [modkit](https://github.com/nanoporetech/modkit)
- [Minimap2](https://github.com/lh3/minimap2)

---

## License
MIT — see [LICENSE](LICENSE)
