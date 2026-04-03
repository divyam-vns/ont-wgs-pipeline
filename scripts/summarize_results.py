#!/usr/bin/env python3
"""
scripts/summarize_results.py
Aggregates QC, alignment, variant calling, and methylation stats
from a completed ONT WGS pipeline run into a single TSV + printed table.

Usage:
    python scripts/summarize_results.py --results results/nextflow --sample DEMO
    python scripts/summarize_results.py --results results/snakemake --sample DEMO
"""

import argparse
import gzip
import os
import re
from pathlib import Path


def parse_flagstat(path: Path) -> dict:
    """Parse samtools flagstat output."""
    stats = {}
    if not path.exists():
        return stats
    with open(path) as f:
        for line in f:
            m = re.match(r"(\d+) \+ (\d+) (.+)", line.strip())
            if m:
                key = m.group(3).split("(")[0].strip().replace(" ", "_")
                stats[key] = int(m.group(1))
    return stats


def parse_nanostats(path: Path) -> dict:
    """Parse NanoPlot NanoStats.txt."""
    stats = {}
    if not path.exists():
        return stats
    with open(path) as f:
        for line in f:
            if ":" in line:
                k, v = line.strip().split(":", 1)
                stats[k.strip()] = v.strip()
    return stats


def count_vcf_variants(vcf_path: Path, sv: bool = False) -> dict:
    """Count PASS variants in a VCF (supports .vcf.gz)."""
    counts = {"TOTAL": 0, "PASS": 0}
    if not vcf_path.exists():
        return counts
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            counts["TOTAL"] += 1
            fields = line.split("\t")
            if len(fields) > 6 and fields[6] in ("PASS", "."):
                counts["PASS"] += 1
            if sv and len(fields) > 7:
                info = fields[7]
                m = re.search(r"SVTYPE=(\w+)", info)
                if m:
                    svtype = m.group(1)
                    counts[svtype] = counts.get(svtype, 0) + 1
    return counts


def count_methylation_sites(bedmethyl_path: Path) -> dict:
    """Count total and highly methylated CpG sites from bedMethyl."""
    stats = {
        "total_CpG_sites": 0,
        "high_meth_sites_gt80pct": 0,
        "low_meth_sites_lt20pct": 0
    }
    if not bedmethyl_path.exists():
        return stats
    opener = gzip.open if str(bedmethyl_path).endswith(".gz") else open
    with opener(bedmethyl_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 11:
                continue
            stats["total_CpG_sites"] += 1
            try:
                fraction = float(fields[10])
                if fraction > 0.80:
                    stats["high_meth_sites_gt80pct"] += 1
                elif fraction < 0.20:
                    stats["low_meth_sites_lt20pct"] += 1
            except (ValueError, IndexError):
                pass
    return stats


def main():
    parser = argparse.ArgumentParser(description="Summarize ONT WGS pipeline results")
    parser.add_argument("--results", required=True, help="Path to results directory")
    parser.add_argument("--sample",  required=True, help="Sample name")
    parser.add_argument("--out",     default=None,   help="Output TSV (default: stdout)")
    args = parser.parse_args()

    R = Path(args.results)
    S = args.sample

    print(f"\n{'='*60}")
    print(f"  ONT WGS Pipeline — Results Summary")
    print(f"  Sample : {S}")
    print(f"  Results: {R}")
    print(f"{'='*60}\n")

    summary = {"sample": S}

    # ---- Alignment ----
    flagstat = parse_flagstat(R / "alignment" / f"{S}.flagstat.txt")
    summary["total_reads"]     = flagstat.get("in_total", "N/A")
    summary["mapped_reads"]    = flagstat.get("mapped", "N/A")
    summary["duplicate_reads"] = flagstat.get("duplicates", "N/A")
    summary["supplementary"]   = flagstat.get("supplementary", "N/A")

    print("ALIGNMENT")
    print(f"  Total reads       : {summary['total_reads']}")
    print(f"  Mapped reads      : {summary['mapped_reads']}")
    print(f"  Supplementary     : {summary['supplementary']}")
    print()

    # ---- NanoPlot QC ----
    nanostats = parse_nanostats(R / "qc" / S / "NanoStats.txt")
    summary["median_read_length"] = nanostats.get("Median read length", "N/A")
    summary["median_qual"]        = nanostats.get("Median read quality", "N/A")
    summary["N50"]                = nanostats.get("Read length N50", "N/A")
    summary["total_bases_Gb"]     = nanostats.get("Total bases", "N/A")

    print("QC (NanoPlot)")
    print(f"  Median read length: {summary['median_read_length']}")
    print(f"  Median quality    : {summary['median_qual']}")
    print(f"  N50               : {summary['N50']}")
    print(f"  Total bases       : {summary['total_bases_Gb']}")
    print()

    # ---- SNP/indel ----
    snp_vcf = R / "variants" / "snp" / S / "merge_output.vcf.gz"
    snp = count_vcf_variants(snp_vcf)
    summary["snp_total"] = snp.get("TOTAL", "N/A")
    summary["snp_pass"]  = snp.get("PASS", "N/A")

    print("VARIANTS — SNP/indel (Clair3)")
    print(f"  Total calls       : {summary['snp_total']}")
    print(f"  PASS calls        : {summary['snp_pass']}")
    print()

    # ---- SV ----
    sv_vcf = R / "variants" / "sv" / f"{S}.sniffles.vcf.gz"
    sv = count_vcf_variants(sv_vcf, sv=True)
    summary["sv_total"] = sv.pop("TOTAL", "N/A")
    summary["sv_pass"]  = sv.pop("PASS", "N/A")
    for svtype, n in sv.items():
        summary[f"sv_{svtype}"] = n

    print("VARIANTS — Structural (Sniffles2)")
    print(f"  Total SV calls    : {summary['sv_total']}")
    print(f"  PASS SV calls     : {summary['sv_pass']}")
    for k, v in summary.items():
        if k.startswith("sv_") and k not in ("sv_total", "sv_pass"):
            print(f"  {k.replace('sv_', ''):18s}: {v}")
    print()

    # ---- Methylation ----
    bed = R / "methylation" / f"{S}.bedmethyl.gz"
    meth = count_methylation_sites(bed)
    summary.update(meth)

    print("METHYLATION (modkit)")
    print(f"  Total CpG sites   : {meth['total_CpG_sites']}")
    print(f"  High meth (>80%)  : {meth['high_meth_sites_gt80pct']}")
    print(f"  Low meth  (<20%)  : {meth['low_meth_sites_lt20pct']}")
    print()

    # ---- Write TSV ----
    out_path = args.out or str(R / "summary" / f"{S}_summary.tsv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    with open(out_path, "w") as f:
        f.write("\t".join(str(k) for k in summary.keys()) + "\n")
        f.write("\t".join(str(v) for v in summary.values()) + "\n")

    print(f"Summary written to: {out_path}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
