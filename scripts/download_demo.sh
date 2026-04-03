#!/usr/bin/env bash
# =============================================================
# scripts/download_demo.sh
# Downloads ONT wf-human-variation demo data + GRCh38 subset
# Source: ONT EPI2ME Labs public S3 bucket (no auth required)
# =============================================================
set -euo pipefail

DEMO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data/demo"
mkdir -p "${DEMO_DIR}"

echo "============================================"
echo " ONT WGS Pipeline — Demo Data Downloader"
echo "============================================"
echo "Destination: ${DEMO_DIR}"
echo ""

# ---- 1. ONT wf-human-variation demo BAM + reference ----
DEMO_URL="https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-human-variation/demo_data.tar.gz"
DEMO_TAR="${DEMO_DIR}/demo_data.tar.gz"

if [ ! -f "${DEMO_DIR}/demo.bam" ]; then
    echo "[1/3] Downloading ONT demo data (~300 MB)..."
    echo "      Source: ${DEMO_URL}"
    wget --progress=bar:force:noscroll \
         -O "${DEMO_TAR}" \
         "${DEMO_URL}"

    echo "[1/3] Extracting demo archive..."
    tar -xzvf "${DEMO_TAR}" -C "${DEMO_DIR}" --strip-components=1
    rm -f "${DEMO_TAR}"
    echo "[1/3] Demo BAM extracted."
else
    echo "[1/3] Demo BAM already present — skipping download."
fi

# ---- 2. Index the demo BAM if needed ----
if [ ! -f "${DEMO_DIR}/demo.bam.bai" ]; then
    echo "[2/3] Indexing demo BAM..."
    samtools index "${DEMO_DIR}/demo.bam"
    echo "[2/3] BAM indexed."
else
    echo "[2/3] BAM index already present — skipping."
fi

# ---- 3. Index the reference FASTA if needed ----
if [ -f "${DEMO_DIR}/demo.fasta" ] && [ ! -f "${DEMO_DIR}/demo.fasta.fai" ]; then
    echo "[3/3] Indexing reference FASTA..."
    samtools faidx "${DEMO_DIR}/demo.fasta"
    echo "[3/3] Reference indexed."
else
    echo "[3/3] Reference index already present — skipping."
fi

# ---- Summary ----
echo ""
echo "============================================"
echo " Demo data ready:"
ls -lh "${DEMO_DIR}/"
echo ""
echo " Chemistry : R10.4.1 e8.2 400bps HAC"
echo " Basecaller: dna_r10.4.1_e8.2_400bps_hac@v4.2.0"
echo " Mod cfg   : dna_r10.4.1_e8.2_400bps_hac@v4.2.0_5mCG_5hmCG@v2"
echo " BAM tags  : MM/ML (5mCG + 5hmCG methylation)"
echo ""
echo " Next steps:"
echo "   Nextflow : nextflow run nextflow/main.nf --input data/demo/demo.bam --ref data/demo/demo.fasta --bed data/demo/demo.bed --sample DEMO --outdir results/nextflow -profile podman"
echo "   Snakemake: snakemake -s snakemake/Snakefile --configfile snakemake/config.yaml --cores 8 --use-conda -n"
echo "============================================"
