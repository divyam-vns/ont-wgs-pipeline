#!/usr/bin/env bash
# =============================================================
# scripts/setup_env.sh
# Sets up all dependencies for the ONT WGS pipeline on macOS
# Stack: Podman (containerized Linux) + Miniconda3 + Nextflow
# =============================================================
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_ENVS=("ont_align" "ont_variants" "ont_methylation")

echo "============================================"
echo " ONT WGS Pipeline — Environment Setup"
echo " Platform: $(uname -s) $(uname -m)"
echo "============================================"
echo ""

# ---- Helper ----
check_cmd() {
    if command -v "$1" &>/dev/null; then
        echo "  [✓] $1 found: $(command -v $1)"
    else
        echo "  [✗] $1 NOT found"
        return 1
    fi
}

# ---- 1. Check Podman ----
echo "[1/5] Checking Podman..."
if ! check_cmd podman; then
    echo ""
    echo "  Podman not found. Install via Homebrew:"
    echo "    brew install podman"
    echo "    podman machine init"
    echo "    podman machine start"
    echo ""
    echo "  Or download from: https://podman.io/getting-started/installation"
    echo ""
    read -rp "  Continue without Podman? (y/n): " yn
    [[ "$yn" == "y" ]] || exit 1
else
    # Check if podman machine is running
    if podman machine list 2>/dev/null | grep -q "Running"; then
        echo "  [✓] Podman machine is running"
    else
        echo "  [!] Podman machine not running. Starting..."
        podman machine start || echo "  [!] Could not start Podman machine — start manually with: podman machine start"
    fi
fi

# ---- 2. Check / install Nextflow ----
echo ""
echo "[2/5] Checking Nextflow..."
if ! check_cmd nextflow; then
    echo "  Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    mkdir -p "$HOME/.local/bin"
    mv nextflow "$HOME/.local/bin/"
    export PATH="$HOME/.local/bin:$PATH"
    echo "  [✓] Nextflow installed at $HOME/.local/bin/nextflow"
    echo "  Add to your shell profile: export PATH=\"\$HOME/.local/bin:\$PATH\""
else
    NF_VER=$(nextflow -version 2>&1 | head -1)
    echo "  Version: ${NF_VER}"
fi

# ---- 3. Check Conda ----
echo ""
echo "[3/5] Checking Conda (Miniconda3)..."
if ! check_cmd conda; then
    echo "  Conda not found."
    echo "  Install Miniconda3 from: https://docs.conda.io/en/latest/miniconda.html"
    echo "  Or via Homebrew: brew install miniconda"
    exit 1
else
    CONDA_VER=$(conda --version)
    echo "  Version: ${CONDA_VER}"
fi

# ---- 4. Create conda environments ----
echo ""
echo "[4/5] Creating Conda environments..."
for ENV in "${CONDA_ENVS[@]}"; do
    ENV_FILE="${REPO_ROOT}/snakemake/envs/${ENV#ont_}.yaml"
    # map env name to yaml: ont_align → align.yaml etc.
    case "${ENV}" in
        ont_align)       YAML="${REPO_ROOT}/snakemake/envs/align.yaml" ;;
        ont_variants)    YAML="${REPO_ROOT}/snakemake/envs/variants.yaml" ;;
        ont_methylation) YAML="${REPO_ROOT}/snakemake/envs/methylation.yaml" ;;
    esac

    if conda env list | grep -q "^${ENV} "; then
        echo "  [✓] ${ENV} already exists — skipping"
    else
        echo "  Creating ${ENV} from ${YAML}..."
        conda env create -f "${YAML}" -n "${ENV}" --quiet
        echo "  [✓] ${ENV} created"
    fi
done

# ---- 5. Check Snakemake ----
echo ""
echo "[5/5] Checking Snakemake..."
if ! check_cmd snakemake; then
    echo "  Installing Snakemake into base conda..."
    conda install -c bioconda -c conda-forge snakemake --quiet -y
fi
check_cmd snakemake
echo "  Version: $(snakemake --version)"

# ---- Pull key Docker/Podman images ----
echo ""
echo "Pulling container images (optional but speeds up first run)..."
IMAGES=(
    "quay.io/biocontainers/minimap2:2.26--he4a0461_2"
    "quay.io/biocontainers/sniffles:2.3.3--pyhdfd78af_0"
    "quay.io/biocontainers/nanoplot:1.42.0--pyhdfd78af_0"
    "docker.io/ontresearch/modkit:latest"
    "docker.io/hkubal/clair3:latest"
)
for IMG in "${IMAGES[@]}"; do
    echo "  Pulling ${IMG}..."
    podman pull "${IMG}" 2>/dev/null && echo "  [✓] ${IMG}" || echo "  [!] Could not pull ${IMG} — will pull on first run"
done

# ---- Done ----
echo ""
echo "============================================"
echo " Setup complete!"
echo ""
echo " Quick start:"
echo "   bash scripts/download_demo.sh"
echo ""
echo "   # Nextflow (Podman):"
echo "   nextflow run nextflow/main.nf \\"
echo "     --input data/demo/demo.bam \\"
echo "     --ref   data/demo/demo.fasta \\"
echo "     --bed   data/demo/demo.bed \\"
echo "     --sample DEMO --outdir results/nextflow -profile podman"
echo ""
echo "   # Snakemake:"
echo "   snakemake -s snakemake/Snakefile \\"
echo "     --configfile snakemake/config.yaml \\"
echo "     --cores 8 --use-conda -n"
echo "============================================"
