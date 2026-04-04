#!/usr/bin/env python3
import argparse
import gzip
import os
import re
from pathlib import Path


def parse_flagstat(path):
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


def parse_nanostats(path):
    stats = {}
    if not path.exists():
        return stats
    with open(path) as f:
        for line in f:
            if ":" in line:
                k, v = line.strip().split(":", 1)
                stats[k.strip()] = v.strip()
    return stats


def count_vcf_variants(vcf_path, sv=False):
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


def count_methylation_sites(bedmethyl_path):
    stats = {"total_CpG_sites": 0, "high_meth_gt80": 0, "low_meth_lt20": 0}
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
                    stats["high_meth_gt80"] += 1
                elif fraction < 0.20:
                    stats["low_meth_lt20"] += 1
            except (ValueError, IndexError):
                pass
    return stats


def main():
    parser = argparse.ArgumentParser(description="Summarize ONT WGS pipeline results")
    parser.add_argument("--results", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    R = Path(args.results)
    S = args.sample

    print("=" * 60)
    print("ONT WGS Pipeline Results Summary")
    print("Sample : " + S)
    print("Results: " + str(R))
    print("=" * 60)

    summary = {"sample": S}

    flagstat = parse_flagstat(R / "alignment" / (S + ".flagstat.txt"))
    summary["total_reads"] = flagstat.get("in_total", "N/A")
    summary["mapped_reads"] = flagstat.get("mapped", "N/A")
    summary["supplementary"] = flagstat.get("supplementary", "N/A")

    print("ALIGNMENT")
    print("  Total reads  : " + str(summary["total_reads"]))
    print("  Mapped reads : " + str(summary["mapped_reads"]))

    nanostats = parse_nanostats(R / "qc" / S / "NanoStats.txt")
    summary["median_read_length"] = nanostats.get("Median read length", "N/A")
    summary["median_qual"] = nanostats.get("Median read quality", "N/A")
    summary["N50"] = nanostats.get("Read length N50", "N/A")
    summary["total_bases"] = nanostats.get("Total bases", "N/A")

    print("QC")
    print("  Median length: " + str(summary["median_read_length"]))
    print("  Median qual  : " + str(summary["median_qual"]))
    print("  N50          : " + str(summary["N50"]))

    snp_vcf = R / "variants" / "snp" / S / "merge_output.vcf.gz"
    snp = count_vcf_variants(snp_vcf)
    summary["snp_total"] = snp.get("TOTAL", "N/A")
    summary["snp_pass"] = snp.get("PASS", "N/A")

    print("SNP/indel (Clair3)")
    print("  Total: " + str(summary["snp_total"]))
    print("  PASS : " + str(summary["snp_pass"]))

    sv_vcf = R / "variants" / "sv" / (S + ".sniffles.vcf.gz")
    sv = count_vcf_variants(sv_vcf, sv=True)
    summary["sv_total"] = sv.pop("TOTAL", "N/A")
    summary["sv_pass"] = sv.pop("PASS", "N/A")
    for svtype, n in sv.items():
        summary["sv_" + svtype] = n

    print("SV (Sniffles2)")
    print("  Total: " + str(summary["sv_total"]))
    print("  PASS : " + str(summary["sv_pass"]))

    bed = R / "methylation" / (S + ".bedmethyl.gz")
    meth = count_methylation_sites(bed)
    summary.update(meth)

    print("Methylation (modkit)")
    print("  CpG sites  : " + str(meth["total_CpG_sites"]))
    print("  High (>80%): " + str(meth["high_meth_gt80"]))
    print("  Low  (<20%): " + str(meth["low_meth_lt20"]))

    out_path = args.out or str(R / "summary" / (S + "_summary.tsv"))
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        f.write("\t".join(str(k) for k in summary.keys()) + "\n")
        f.write("\t".join(str(v) for v in summary.values()) + "\n")
    print("Written to: " + out_path)


if __name__ == "__main__":
    main()
