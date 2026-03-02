#!/usr/bin/env python3
"""
summarize_variants.py
─────────────────────
Parse annotated VCF (SnpEff output) and produce:
  - TSV summary table
  - PDF report with statistics and plots
"""

import sys
import os
import re
import gzip
from collections import defaultdict, Counter
import json

def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def parse_info(info_str):
    """Parse VCF INFO field into dict."""
    d = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            d[k] = v
        else:
            d[field] = True
    return d

def parse_ann(ann_str):
    """Parse SnpEff ANN field — returns list of effect dicts."""
    effects = []
    for allele in ann_str.split(","):
        parts = allele.split("|")
        if len(parts) >= 4:
            effects.append({
                "allele":   parts[0],
                "effect":   parts[1],
                "impact":   parts[2],
                "gene":     parts[3],
                "feature":  parts[4] if len(parts) > 4 else "",
                "hgvs_c":   parts[9] if len(parts) > 9 else "",
                "hgvs_p":   parts[10] if len(parts) > 10 else "",
            })
    return effects

def classify_variant(ref, alt):
    """Classify as SNP, MNP, INS, DEL."""
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    elif len(ref) == len(alt):
        return "MNP"
    elif len(ref) < len(alt):
        return "INS"
    else:
        return "DEL"

def ts_tv(ref, alt):
    """Transition or Transversion."""
    transitions = {("A","G"),("G","A"),("C","T"),("T","C")}
    pair = (ref.upper(), alt.upper())
    if pair in transitions:
        return "Ti"
    return "Tv"

def main():
    if len(sys.argv) < 3:
        print("Usage: summarize_variants.py <input.vcf[.gz]> <output_prefix>")
        sys.exit(1)

    vcf_in  = sys.argv[1]
    out_pfx = sys.argv[2]

    variants = []
    samples  = []
    impacts  = Counter()
    effects  = Counter()
    var_types = Counter()
    genes    = Counter()
    ti_count = 0
    tv_count = 0

    print(f"📂 Parsing VCF: {vcf_in}")

    with open_vcf(vcf_in) as fh:
        for line in fh:
            line = line.rstrip("\n")

            # Header
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.lstrip("#").split("\t")
                if len(cols) > 9:
                    samples = cols[9:]
                continue

            # Variant line
            parts = line.split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, vid, ref, alt_field, qual, filt, info_str = parts[:8]
            fmt  = parts[8] if len(parts) > 8 else ""
            gts  = parts[9:] if len(parts) > 9 else []

            # Skip filtered
            if filt not in ("PASS", "."):
                continue

            info = parse_info(info_str)
            alts = alt_field.split(",")

            for alt in alts:
                vtype = classify_variant(ref, alt)
                var_types[vtype] += 1

                # Ti/Tv for SNPs
                if vtype == "SNP":
                    label = ts_tv(ref, alt)
                    if label == "Ti":
                        ti_count += 1
                    else:
                        tv_count += 1

                # SnpEff annotations
                ann_effects = []
                if "ANN" in info:
                    ann_effects = parse_ann(info["ANN"])
                    for eff in ann_effects[:1]:  # primary effect
                        impacts[eff["impact"]] += 1
                        effects[eff["effect"]] += 1
                        if eff["gene"]:
                            genes[eff["gene"]] += 1

                # AC/AF from INFO
                ac = info.get("AC", ".")
                af = info.get("AF", ".")

                row = {
                    "CHROM": chrom, "POS": pos, "ID": vid,
                    "REF": ref, "ALT": alt, "QUAL": qual,
                    "FILTER": filt, "TYPE": vtype,
                    "AC": ac, "AF": af,
                    "IMPACT": ann_effects[0]["impact"] if ann_effects else ".",
                    "EFFECT": ann_effects[0]["effect"] if ann_effects else ".",
                    "GENE":   ann_effects[0]["gene"]   if ann_effects else ".",
                    "HGVS_C": ann_effects[0]["hgvs_c"] if ann_effects else ".",
                    "HGVS_P": ann_effects[0]["hgvs_p"] if ann_effects else ".",
                }
                variants.append(row)

    n_total = len(variants)
    ti_tv_ratio = round(ti_count / tv_count, 3) if tv_count > 0 else "N/A"

    print(f"  Total PASS variants : {n_total:,}")
    print(f"  SNPs                : {var_types['SNP']:,}")
    print(f"  INDELs (INS+DEL)    : {var_types['INS'] + var_types['DEL']:,}")
    print(f"  Ti/Tv ratio         : {ti_tv_ratio}")
    print(f"  HIGH impact         : {impacts.get('HIGH', 0):,}")
    print(f"  MODERATE impact     : {impacts.get('MODERATE', 0):,}")

    # ── Write TSV ──────────────────────────────────────────────────────────
    tsv_path = f"{out_pfx}_variants.tsv"
    if variants:
        headers = list(variants[0].keys())
        with open(tsv_path, "w") as f:
            f.write("\t".join(headers) + "\n")
            for v in variants:
                f.write("\t".join(str(v[h]) for h in headers) + "\n")
        print(f"  ✓ TSV written: {tsv_path}")

    # ── Write JSON stats ───────────────────────────────────────────────────
    stats = {
        "total_variants": n_total,
        "by_type": dict(var_types),
        "ti_count": ti_count,
        "tv_count": tv_count,
        "ti_tv_ratio": ti_tv_ratio,
        "by_impact": dict(impacts),
        "top_effects": dict(effects.most_common(10)),
        "top_genes": dict(genes.most_common(20)),
        "samples": samples,
    }
    json_path = f"{out_pfx}_stats.json"
    with open(json_path, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"  ✓ Stats JSON: {json_path}")

    # ── Generate PDF report ────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.gridspec import GridSpec
        from reportlab.lib.pagesizes import A4
        from reportlab.lib import colors
        from reportlab.lib.units import cm
        from reportlab.platypus import (SimpleDocTemplate, Table, TableStyle,
                                         Paragraph, Spacer, Image, HRFlowable)
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.enums import TA_CENTER

        pdf_path = f"{out_pfx}_report.pdf"
        fig_dir  = os.path.dirname(out_pfx) or "."

        # --- Figure 1: Variant types pie ---
        fig1, axes = plt.subplots(1, 3, figsize=(14, 4.5))
        fig1.patch.set_facecolor("#f8f9fa")

        # Pie: variant types
        labels1 = [k for k, v in var_types.items() if v > 0]
        vals1   = [var_types[k] for k in labels1]
        clrs1   = ["#4e79a7", "#f28e2b", "#e15759", "#76b7b2"]
        axes[0].pie(vals1, labels=labels1, colors=clrs1[:len(labels1)],
                    autopct="%1.1f%%", startangle=90,
                    wedgeprops={"edgecolor": "white", "linewidth": 1.5})
        axes[0].set_title("Variant Types", fontweight="bold", fontsize=12)

        # Bar: impact
        impact_order = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
        imp_vals = [impacts.get(i, 0) for i in impact_order]
        imp_clrs = ["#e15759", "#f28e2b", "#59a14f", "#76b7b2"]
        bars = axes[1].bar(impact_order, imp_vals, color=imp_clrs, edgecolor="white")
        axes[1].set_title("Functional Impact", fontweight="bold", fontsize=12)
        axes[1].set_ylabel("Count")
        for bar, val in zip(bars, imp_vals):
            if val > 0:
                axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                             str(val), ha="center", va="bottom", fontsize=9)

        # Bar: top effects
        top_eff = effects.most_common(8)
        if top_eff:
            eff_names = [e[0].replace("_", " ")[:20] for e, _ in top_eff]
            eff_vals  = [v for _, v in top_eff]
            axes[2].barh(eff_names[::-1], eff_vals[::-1], color="#4e79a7", edgecolor="white")
            axes[2].set_title("Top Variant Effects", fontweight="bold", fontsize=12)
            axes[2].set_xlabel("Count")

        plt.tight_layout(pad=2)
        fig1_path = os.path.join(fig_dir, "fig_variants.png")
        fig1.savefig(fig1_path, dpi=150, bbox_inches="tight")
        plt.close(fig1)

        # --- Figure 2: Top genes ---
        fig2, ax = plt.subplots(figsize=(10, 4))
        fig2.patch.set_facecolor("#f8f9fa")
        top_g = genes.most_common(15)
        if top_g:
            g_names = [g for g, _ in top_g]
            g_vals  = [v for _, v in top_g]
            ax.bar(g_names, g_vals, color="#4e79a7", edgecolor="white")
            ax.set_title("Most Mutated Genes", fontweight="bold", fontsize=12)
            ax.set_ylabel("Variant Count")
            ax.tick_params(axis="x", rotation=45)
        plt.tight_layout()
        fig2_path = os.path.join(fig_dir, "fig_genes.png")
        fig2.savefig(fig2_path, dpi=150, bbox_inches="tight")
        plt.close(fig2)

        # --- Build PDF ---
        doc  = SimpleDocTemplate(pdf_path, pagesize=A4,
                                  topMargin=1.5*cm, bottomMargin=1.5*cm,
                                  leftMargin=2*cm, rightMargin=2*cm)
        styles = getSampleStyleSheet()
        title_style = ParagraphStyle("title", parent=styles["Title"],
                                      fontSize=18, spaceAfter=6,
                                      textColor=colors.HexColor("#2c3e50"))
        h1 = ParagraphStyle("h1", parent=styles["Heading1"],
                             fontSize=13, textColor=colors.HexColor("#2980b9"),
                             spaceBefore=12, spaceAfter=4)
        normal = styles["Normal"]

        story = []
        story.append(Paragraph("🧬 Variant Calling Pipeline — Summary Report", title_style))
        story.append(Paragraph("Germline SNP & INDEL Analysis | GATK4 Best Practices", normal))
        story.append(HRFlowable(width="100%", thickness=1,
                                 color=colors.HexColor("#2980b9"), spaceAfter=12))

        # Key stats table
        story.append(Paragraph("Key Statistics", h1))
        kv = [
            ["Metric", "Value"],
            ["Total PASS Variants", f"{n_total:,}"],
            ["SNPs", f"{var_types.get('SNP',0):,}"],
            ["Insertions", f"{var_types.get('INS',0):,}"],
            ["Deletions", f"{var_types.get('DEL',0):,}"],
            ["Ti/Tv Ratio", str(ti_tv_ratio)],
            ["HIGH Impact Variants", f"{impacts.get('HIGH',0):,}"],
            ["MODERATE Impact Variants", f"{impacts.get('MODERATE',0):,}"],
            ["Samples Analyzed", str(len(samples))],
        ]
        tbl = Table(kv, colWidths=[9*cm, 6*cm])
        tbl.setStyle(TableStyle([
            ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#2980b9")),
            ("TEXTCOLOR",  (0,0), (-1,0), colors.white),
            ("FONTNAME",   (0,0), (-1,0), "Helvetica-Bold"),
            ("FONTSIZE",   (0,0), (-1,-1), 10),
            ("ROWBACKGROUNDS", (0,1), (-1,-1),
             [colors.HexColor("#ecf0f1"), colors.white]),
            ("GRID", (0,0), (-1,-1), 0.5, colors.grey),
            ("ALIGN", (1,0), (1,-1), "CENTER"),
            ("TOPPADDING",    (0,0), (-1,-1), 5),
            ("BOTTOMPADDING", (0,0), (-1,-1), 5),
        ]))
        story.append(tbl)
        story.append(Spacer(1, 0.4*cm))

        # Figures
        story.append(Paragraph("Variant Distribution & Functional Impact", h1))
        story.append(Image(fig1_path, width=16*cm, height=5.5*cm))
        story.append(Spacer(1, 0.3*cm))

        if os.path.exists(fig2_path) and genes:
            story.append(Paragraph("Most Mutated Genes", h1))
            story.append(Image(fig2_path, width=16*cm, height=4.5*cm))

        # Pipeline info
        story.append(Spacer(1, 0.4*cm))
        story.append(Paragraph("Pipeline Information", h1))
        pipe_data = [
            ["Step", "Tool", "Version"],
            ["Quality Control",    "FastQC",               "0.12"],
            ["Read Trimming",      "fastp",                "0.23"],
            ["Alignment",          "BWA-MEM2",             "2.2.1"],
            ["Duplicate Marking",  "GATK MarkDuplicates",  "4.4"],
            ["Base Recalibration", "GATK BQSR",            "4.4"],
            ["Variant Calling",    "GATK HaplotypeCaller", "4.4"],
            ["Joint Genotyping",   "GATK GenotypeGVCFs",   "4.4"],
            ["Hard Filtering",     "GATK VariantFiltration","4.4"],
            ["Annotation",         "SnpEff",               "5.2"],
            ["QC Aggregation",     "MultiQC",              "1.19"],
        ]
        tbl2 = Table(pipe_data, colWidths=[6*cm, 6*cm, 4*cm])
        tbl2.setStyle(TableStyle([
            ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#27ae60")),
            ("TEXTCOLOR",  (0,0), (-1,0), colors.white),
            ("FONTNAME",   (0,0), (-1,0), "Helvetica-Bold"),
            ("FONTSIZE",   (0,0), (-1,-1), 9),
            ("ROWBACKGROUNDS", (0,1), (-1,-1),
             [colors.HexColor("#eafaf1"), colors.white]),
            ("GRID", (0,0), (-1,-1), 0.5, colors.grey),
            ("TOPPADDING",    (0,0), (-1,-1), 4),
            ("BOTTOMPADDING", (0,0), (-1,-1), 4),
        ]))
        story.append(tbl2)

        doc.build(story)
        print(f"  ✓ PDF report: {pdf_path}")

    except ImportError as e:
        print(f"  ⚠ PDF skipped (missing library: {e}). TSV+JSON still written.")

    print("\n✅ Variant summary complete.")

if __name__ == "__main__":
    main()
