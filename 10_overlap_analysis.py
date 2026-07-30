#!/usr/bin/env python
"""
10_overlap_analysis.py

This script aims to connect Memory B and Plasma cell findings:
1. TF overlap from pySCENIC
2. Shared vs unique Lasso predictors per TF with direction comparison
3. DE gene overlap with Lasso predictors
4. pySCENIC target overlap for shared TFs
5. TF recurrence heatmap
6. g:Profiler enrichment on shared vs unique gene sets

Usage:
    python 10_overlap_analysis.py \
        --pyscenic_dir  /path/to/pyscenic_out \
        --lasso_dir     /path/to/lasso_out \
        --de_csv        /path/to/pseudobulk_out/de_memoryB_vs_plasma.csv \
        --pairs_csv     /path/to/grn_analysis_out/most_common_pairs.csv \
        --outdir        /path/to/overlap_out \
        --tfs           BACH2 ATF4 XBP1
"""
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from gprofiler import GProfiler


def clean_genes(gene_set):
    """Remove NaN/float values from gene sets."""
    return {g for g in gene_set
            if isinstance(g, str) and g.strip() not in ("", "nan")}


def load_lasso_coefs(lasso_dir, tf, cell_type):
    path = f"{lasso_dir}/{tf}/{cell_type}_{tf}_lasso_coefs.csv"
    if not os.path.exists(path):
        return pd.Series(dtype=float)
    df = pd.read_csv(path, index_col=0, header=None)
    df.columns = ["coefficient"]
    df["coefficient"] = pd.to_numeric(df["coefficient"], errors="coerce")
    df = df.dropna()
    return df["coefficient"]


def load_recurrent_targets(pairs_csv, cell_type, min_donors=6):
    pairs = pd.read_csv(pairs_csv)
    if "pair" in pairs.columns:
        pairs[["TF", "target"]] = pairs["pair"].str.split("->", expand=True)
    sub = pairs[pairs["cell_type"] == cell_type].copy()
    return sub[sub["n_donors"] >= min_donors]


def run_gprofiler(gene_list, label, outdir):
    gene_list = [str(g).strip() for g in gene_list
                 if str(g).strip() not in ("", "nan")]
    if len(gene_list) < 5:
        print(f"  Skipping {label} — too few genes ({len(gene_list)})", flush=True)
        return None
    print(f"  g:Profiler: {label} ({len(gene_list)} genes)", flush=True)
    gp = GProfiler(return_dataframe=True)
    results = gp.profile(
        organism="hsapiens",
        query=gene_list,
        sources=["GO:BP", "GO:MF", "KEGG", "REAC", "WP"],
        significance_threshold_method="fdr",
        user_threshold=0.05
    )
    if results.empty:
        print(f"  No significant terms for {label}", flush=True)
        return results
    safe = label.replace(" ", "_").replace("/", "_")
    results.to_csv(f"{outdir}/{safe}_enrichment.csv", index=False)
    return results


def plot_venn_bar(sets, title, outpath):
    categories = ["Memory B only", "Shared", "Plasma only"]
    values = [len(sets[0] - sets[1]), len(sets[0] & sets[1]), len(sets[1] - sets[0])]
    colours = ["#1565C0", "#9C27B0", "#C62828"]
    fig, ax = plt.subplots(figsize=(8, 4))
    bars = ax.barh(categories, values, color=colours, edgecolor="white", height=0.5)
    for bar, val in zip(bars, values):
        ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
                str(val), va="center", fontsize=12, fontweight="bold")
    ax.set_xlabel("Number of genes / TFs")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlim(0, max(values) * 1.2 if max(values) > 0 else 10)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_direction_scatter(shared_genes, mb_coefs, pl_coefs, tf, outdir):
    shared_genes = clean_genes(shared_genes)
    if len(shared_genes) == 0:
        return
    df = pd.DataFrame({
        "Memory_B": mb_coefs.reindex(list(shared_genes)).fillna(0),
        "Plasma":   pl_coefs.reindex(list(shared_genes)).fillna(0),
    })
    colours = []
    for _, row in df.iterrows():
        if row.Memory_B > 0 and row.Plasma > 0:
            colours.append("#1565C0")
        elif row.Memory_B < 0 and row.Plasma < 0:
            colours.append("#C62828")
        else:
            colours.append("#E91E8C")
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(df.Memory_B, df.Plasma, c=colours, alpha=0.8, s=60, edgecolors="white")
    ax.axhline(0, color="grey", linewidth=0.8, linestyle="--")
    ax.axvline(0, color="grey", linewidth=0.8, linestyle="--")
    top = df.assign(mag=df.Memory_B.abs() + df.Plasma.abs()).nlargest(10, "mag")
    for gene, row in top.iterrows():
        ax.annotate(gene, (row.Memory_B, row.Plasma),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")
    patches = [
        mpatches.Patch(color="#1565C0", label="Both positive"),
        mpatches.Patch(color="#C62828", label="Both negative"),
        mpatches.Patch(color="#E91E8C", label="Discordant"),
    ]
    ax.legend(handles=patches, fontsize=9)
    ax.set_xlabel("Lasso coefficient — Memory B cells")
    ax.set_ylabel("Lasso coefficient — Plasma cells")
    ax.set_title(f"{tf}: shared predictor genes across cell types\n({len(shared_genes)} shared genes)",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(f"{outdir}/{tf}_shared_predictor_directions.png",
                dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main(pyscenic_dir, lasso_dir, de_csv, pairs_csv, outdir, tfs):
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(f"{outdir}/enrichment", exist_ok=True)

    #pySCENIC TF overlap
    print("\n=== 1. pySCENIC TF overlap ===", flush=True)
    mb_tfs = set(load_recurrent_targets(pairs_csv, "Memory_B")["TF"].unique())
    pl_tfs = set(load_recurrent_targets(pairs_csv, "Plasma")["TF"].unique())
    shared_tfs  = mb_tfs & pl_tfs
    mb_only_tfs = mb_tfs - pl_tfs
    pl_only_tfs = pl_tfs - mb_tfs

    print(f"  Memory B TFs: {len(mb_tfs)}, Plasma TFs: {len(pl_tfs)}, Shared: {len(shared_tfs)}", flush=True)
    print(f"  Shared TFs: {sorted(shared_tfs)}", flush=True)
    print(f"  Memory B only: {sorted(mb_only_tfs)}", flush=True)
    print(f"  Plasma only:   {sorted(pl_only_tfs)}", flush=True)

    pd.DataFrame({
        "TF": sorted(mb_tfs | pl_tfs),
        "in_Memory_B": [t in mb_tfs for t in sorted(mb_tfs | pl_tfs)],
        "in_Plasma":   [t in pl_tfs for t in sorted(mb_tfs | pl_tfs)],
        "shared":      [t in shared_tfs for t in sorted(mb_tfs | pl_tfs)]
    }).to_csv(f"{outdir}/pyscenic_tf_overlap.csv", index=False)

    plot_venn_bar([mb_tfs, pl_tfs],
                  "Recurrent TFs: Memory B vs Plasma (pySCENIC)",
                  f"{outdir}/tf_overlap_bar.png")

    #Lasso predictor overlap
    print("\n=== 2. Lasso predictor overlap ===", flush=True)
    lasso_rows = []
    for tf in tfs:
        mb_coefs = load_lasso_coefs(lasso_dir, tf, "Memory_B")
        pl_coefs = load_lasso_coefs(lasso_dir, tf, "Plasma")
        mb_genes = clean_genes(set(mb_coefs.index))
        pl_genes = clean_genes(set(pl_coefs.index))
        shared   = mb_genes & pl_genes
        mb_only  = mb_genes - pl_genes
        pl_only  = pl_genes - mb_genes

        print(f"  {tf}: {len(mb_genes)} MB, {len(pl_genes)} Plasma, {len(shared)} shared", flush=True)

        lasso_rows.append({
            "TF": tf,
            "n_Memory_B": len(mb_genes), "n_Plasma": len(pl_genes),
            "n_shared": len(shared), "n_Memory_B_only": len(mb_only),
            "n_Plasma_only": len(pl_only),
            "shared_genes": ";".join(sorted(shared)),
            "Memory_B_only_genes": ";".join(sorted(mb_only)),
            "Plasma_only_genes": ";".join(sorted(pl_only))
        })

        plot_direction_scatter(shared, mb_coefs, pl_coefs, tf, outdir)
        plot_venn_bar([mb_genes, pl_genes],
                      f"{tf} Lasso predictors: Memory B vs Plasma",
                      f"{outdir}/{tf}_predictor_overlap_bar.png")

        for gene_set, label in [
            (shared,  f"{tf}_shared_predictors"),
            (mb_only, f"{tf}_Memory_B_only_predictors"),
            (pl_only, f"{tf}_Plasma_only_predictors"),
        ]:
            run_gprofiler(list(gene_set), label, f"{outdir}/enrichment")

    pd.DataFrame(lasso_rows).to_csv(
        f"{outdir}/lasso_predictor_overlap_summary.csv", index=False)

    #DE overlap with Lasso predictors 
    print("\n=== 3. DE x Lasso overlap ===", flush=True)
    de = pd.read_csv(de_csv, index_col=0).dropna(subset=["padj"])
    de_sig = de[de["padj"] < 0.05]
    de_up_mb = set(de_sig[de_sig["log2FoldChange"] > 1].index)
    de_up_pl = set(de_sig[de_sig["log2FoldChange"] < -1].index)

    de_overlap_rows = []
    for tf in tfs:
        mb_coefs = load_lasso_coefs(lasso_dir, tf, "Memory_B")
        pl_coefs = load_lasso_coefs(lasso_dir, tf, "Plasma")
        mb_lasso_de = clean_genes(set(mb_coefs.index)) & de_up_mb
        pl_lasso_de = clean_genes(set(pl_coefs.index)) & de_up_pl
        print(f"  {tf}: {len(mb_lasso_de)} MB Lasso+DE, {len(pl_lasso_de)} Plasma Lasso+DE", flush=True)
        if mb_lasso_de:
            print(f"    MB genes: {sorted(mb_lasso_de)}", flush=True)
        if pl_lasso_de:
            print(f"    Plasma genes: {sorted(pl_lasso_de)}", flush=True)
        de_overlap_rows.append({
            "TF": tf,
            "MB_lasso_and_DE_upMB": ";".join(sorted(mb_lasso_de)),
            "n_MB_lasso_and_DE_upMB": len(mb_lasso_de),
            "PL_lasso_and_DE_upPL": ";".join(sorted(pl_lasso_de)),
            "n_PL_lasso_and_DE_upPL": len(pl_lasso_de),
        })

    pd.DataFrame(de_overlap_rows).to_csv(f"{outdir}/lasso_de_overlap.csv", index=False)

    #Shared TF target overlap from pySCENIC
    print("\n=== 4. Shared TF pySCENIC target overlap ===", flush=True)
    if shared_tfs:
        target_rows = []
        for tf in sorted(shared_tfs):
            mb_t = set(load_recurrent_targets(pairs_csv, "Memory_B").query("TF == @tf")["target"])
            pl_t = set(load_recurrent_targets(pairs_csv, "Plasma").query("TF == @tf")["target"])
            shared_t = clean_genes(mb_t) & clean_genes(pl_t)
            print(f"  {tf}: {len(mb_t)} MB, {len(pl_t)} Plasma, {len(shared_t)} shared targets", flush=True)
            target_rows.append({
                "TF": tf, "n_MB": len(mb_t), "n_Plasma": len(pl_t),
                "n_shared": len(shared_t),
                "shared_targets": ";".join(sorted(shared_t))
            })
            if shared_t:
                run_gprofiler(list(shared_t), f"{tf}_shared_pySCENIC_targets",
                              f"{outdir}/enrichment")
        pd.DataFrame(target_rows).to_csv(
            f"{outdir}/pyscenic_shared_tf_targets.csv", index=False)
    else:
        print("  No shared TFs found — skipping target overlap", flush=True)

    #TF recurrence heatmap 
    print("\n=== 5. TF recurrence heatmap ===", flush=True)
    pairs = pd.read_csv(pairs_csv)
    if "pair" in pairs.columns:
        pairs[["TF", "target"]] = pairs["pair"].str.split("->", expand=True)
    all_tfs = sorted(pairs["TF"].unique())
    hmap = pd.DataFrame(index=all_tfs, columns=["Memory_B", "Plasma"], dtype=float).fillna(0)
    for ct in ["Memory_B", "Plasma"]:
        ct_data = pairs[pairs["cell_type"] == ct].groupby("TF")["n_donors"].max()
        hmap[ct] = ct_data
    hmap = hmap.fillna(0)
    hmap = hmap[hmap.sum(axis=1) > 0].sort_values("Memory_B", ascending=False)

    fig, ax = plt.subplots(figsize=(5, max(6, 0.35 * len(hmap))))
    sns.heatmap(hmap, ax=ax, cmap="RdPu", linewidths=0.3, linecolor="white",
                cbar_kws={"label": "Max n donors (recurrence)"},
                annot=True, fmt=".0f", annot_kws={"size": 8})
    ax.set_title("TF recurrence: Memory B vs Plasma\n(max donors with recurrent TF-target pair)",
                 fontsize=10, fontweight="bold")
    ax.set_xticklabels(["Memory B", "Plasma"], rotation=0)
    fig.tight_layout()
    fig.savefig(f"{outdir}/tf_recurrence_heatmap.png",
                dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"\n=== Complete. Outputs in {outdir} ===")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--pyscenic_dir", required=True)
    p.add_argument("--lasso_dir",    required=True)
    p.add_argument("--de_csv",       required=True)
    p.add_argument("--pairs_csv",    required=True)
    p.add_argument("--outdir",       required=True)
    p.add_argument("--tfs", nargs="+", default=["BACH2", "ATF4", "XBP1"])
    args = p.parse_args()
    main(args.pyscenic_dir, args.lasso_dir, args.de_csv,
         args.pairs_csv, args.outdir, args.tfs)
