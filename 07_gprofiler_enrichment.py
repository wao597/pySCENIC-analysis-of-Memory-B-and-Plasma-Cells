#!/usr/bin/env python
"""
07_gprofiler_enrichment.py

This script runs g:Profiler gene set enrichment on four gene lists from the pipeline:
  1. Non-zero Lasso coefficients per TF per cell type
  2. Top DE genes from pseudobulk analysis (up in Memory B, up in Plasma)
  3. Target genes of most recurrent pySCENIC regulons per cell type

For each gene list, queries g:Profiler (GO:BP, GO:MF, KEGG, Reactome,
WikiPathways) and writes a results CSV per gene list & a dot plot of top enriched terms.

Usage:
    python 07_gprofiler_enrichment.py \
        --lasso_dir /path/to/lasso_out \
        --de_csv /path/to/pseudobulk_out/de_memoryB_vs_plasma.csv \
        --grn_csv /path/to/grn_analysis_out/most_common_pairs.csv \
        --outdir /path/to/gprofiler_out \
        --tfs BACH2 ATF4 XBP1 RORA\
        --n_de 200 \
        --top_terms 15
"""


import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from gprofiler import GProfiler


def clean_gene_list(gene_list):
    """Clean gene list by removing NaN, empty strings, and non-string values."""
    cleaned = []
    for g in gene_list:
        # Skip if NaN or None
        if pd.isna(g) or g is None:
            continue
        # Convert to string and strip
        g_str = str(g).strip()
        # Skip empty strings or 'nan' strings
        if g_str == "" or g_str.lower() == "nan" or g_str == "None":
            continue
        # Skip if it starts with a number or contains invalid characters
        
        if g_str[0].isdigit():
            continue
        cleaned.append(g_str)
    return cleaned


def run_gprofiler(gene_list, label, outdir, top_terms=15):
    """Run g:Profiler on a gene list and save results + dot plot."""
    gene_list = clean_gene_list(gene_list)
    
    if len(gene_list) == 0:
        print(f"  WARNING: empty gene list for {label}, skipping", flush=True)
        return None
    
    # Remove duplicates 
    gene_list = list(dict.fromkeys(gene_list))
    
    print(f"  Running g:Profiler: {label} ({len(gene_list)} genes)", flush=True)

    gp = GProfiler(return_dataframe=True)
    
    try:
        results = gp.profile(
            organism="hsapiens",
            query=gene_list,
            sources=["GO:BP", "GO:MF", "KEGG", "REAC", "WP"],
            significance_threshold_method="fdr",
            user_threshold=0.05,
            no_evidences=False
        )
    except Exception as e:
        print(f"  ERROR running g:Profiler for {label}: {e}", flush=True)
        return None

    safe_label = label.replace(" ", "_").replace("/", "_").replace("-", "_")

    if results is None or results.empty:
        print(f"  No significant terms found for {label}", flush=True)
        return results

    results.to_csv(f"{outdir}/{safe_label}_gprofiler_results.csv", index=False)
    print(f"  {len(results)} significant terms found", flush=True)

    # dot plot
    plot_df = results.sort_values("p_value").head(top_terms).copy()
    plot_df["gene_ratio"] = plot_df["intersection_size"] / plot_df["query_size"]
    plot_df["-log10_padj"] = -np.log10(plot_df["p_value"].clip(lower=1e-300))
    plot_df["name"] = plot_df["name"].str[:60]

    fig, ax = plt.subplots(figsize=(10, max(6, 0.5 * len(plot_df))))
    
    scatter = ax.scatter(
        plot_df["gene_ratio"],
        range(len(plot_df)),
        c=plot_df["-log10_padj"],      
        s=plot_df["intersection_size"] * 15,  
        cmap="RdPu",                   
        vmin=0,                        
        vmax=plot_df["-log10_padj"].max(),
        alpha=0.92,
        edgecolors="grey",
        linewidths=0.4
    )
    
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["name"], fontsize=10)
    ax.set_xlabel("Gene ratio (intersection / query size)", fontsize=11)
    ax.set_title(f"Top enriched terms: {label}", fontsize=12, fontweight="bold")
    
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.4, pad=0.02)
    cbar.set_label("-log10(adjusted p-value)\nDarker = more significant",
                   fontsize=10)
    cbar.ax.tick_params(labelsize=9)
    
    for size in [5, 20, 50]:
        ax.scatter([], [], s=size * 15, c="grey", alpha=0.6, label=f"n={size} genes")
    ax.legend(title="Intersection size", loc="lower right", fontsize=8)
    
    fig.tight_layout()ize=7)
    fig.tight_layout()
    fig.savefig(f"{outdir}/{safe_label}_dotplot.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    return results


def load_lasso_genes(lasso_dir, tf, cell_type, direction=None):
    coef_path = f"{lasso_dir}/{tf}/{cell_type}_{tf}_lasso_coefs.csv"
    if not os.path.exists(coef_path):
        print(f"  WARNING: not found: {coef_path}", flush=True)
        return []
    
    try:
        coefs = pd.read_csv(coef_path, index_col=0, header=None)
        coefs.columns = ["coefficient"]
        coefs["coefficient"] = pd.to_numeric(coefs["coefficient"], errors="coerce")
        coefs = coefs.dropna(subset=["coefficient"])
        
        if direction == "positive":
            coefs = coefs[coefs["coefficient"] > 0]
        elif direction == "negative":
            coefs = coefs[coefs["coefficient"] < 0]
        
        genes = [str(g).strip() for g in coefs.index.tolist()]
        return genes
    except Exception as e:
        print(f"  ERROR loading {coef_path}: {e}", flush=True)
        return []


def load_de_genes(de_csv, n=200, direction="up_memory"):
    try:
        de = pd.read_csv(de_csv, index_col=0)
        de = de.dropna(subset=["padj"])
        de = de[de["padj"] < 0.05]
        
        if direction == "up_memory":
            de = de[de["log2FoldChange"] > 0].sort_values("log2FoldChange", ascending=False)
        else:  # up_plasma
            de = de[de["log2FoldChange"] < 0].sort_values("log2FoldChange", ascending=True)
        
        genes = de.head(n).index.tolist()
        return genes
    except Exception as e:
        print(f"  ERROR loading DE genes: {e}", flush=True)
        return []


def load_regulon_targets(grn_csv, tf, cell_type):
    try:
        pairs = pd.read_csv(grn_csv)
        sub = pairs[(pairs["cell_type"] == cell_type) &
                    (pairs["pair"].str.startswith(f"{tf}->"))]
        targets = sub["pair"].str.replace(f"{tf}->", "", regex=False).tolist()
        return targets
    except Exception as e:
        print(f"  ERROR loading regulon targets: {e}", flush=True)
        return []


def main(lasso_dir, de_csv, grn_csv, outdir, tfs, n_de, top_terms):
    os.makedirs(outdir, exist_ok=True)
    all_results = {}

    #Lasso enrichment
    print("\n=== 1. Lasso gene list enrichment ===", flush=True)
    for tf in tfs:
        for cell_type in ["Memory_B", "Plasma"]:
            for direction in [None, "positive", "negative"]:
                genes = load_lasso_genes(lasso_dir, tf, cell_type, direction)
                suffix = f"_{direction}" if direction else "_all"
                label = f"Lasso_{tf}_{cell_type}{suffix}"
                r = run_gprofiler(genes, label, outdir, top_terms)
                if r is not None:
                    all_results[label] = r

    #DE enrichment
    print("\n=== 2. Pseudobulk DE enrichment ===", flush=True)
    if os.path.exists(de_csv):
        for direction, label in [
            ("up_memory", f"DE_top{n_de}_up_in_MemoryB"),
            ("up_plasma",  f"DE_top{n_de}_up_in_Plasma")
        ]:
            genes = load_de_genes(de_csv, n=n_de, direction=direction)
            r = run_gprofiler(genes, label, outdir, top_terms)
            if r is not None:
                all_results[label] = r
    else:
        print(f"  WARNING: DE CSV not found: {de_csv}", flush=True)

    #Regulon target enrichment
    print("\n=== 3. pySCENIC regulon target enrichment ===", flush=True)
    if os.path.exists(grn_csv):
        for tf, cell_type in [("BACH2", "Memory_B"), ("ATF4", "Plasma"), ("XBP1", "Plasma")]:
            genes = load_regulon_targets(grn_csv, tf, cell_type)
            if genes:
                label = f"Regulon_{tf}_{cell_type}_targets"
                r = run_gprofiler(genes, label, outdir, top_terms)
                if r is not None:
                    all_results[label] = r
            else:
                print(f"  No regulon targets found for {tf} in {cell_type}", flush=True)
    else:
        print(f"  WARNING: GRN CSV not found: {grn_csv}", flush=True)

    # Writing a summary table - top 5 terms per analysis
    summary_rows = []
    for label, res in all_results.items():
        if res is not None and not res.empty:
            for _, row in res.sort_values("p_value").head(5).iterrows():
                summary_rows.append({
                    "analysis": label,
                    "source": row["source"],
                    "term_name": row["name"],
                    "padj": row["p_value"],
                    "intersection_size": row["intersection_size"],
                    "gene_ratio": round(row["intersection_size"] / row["query_size"], 3)
                })
    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(
            f"{outdir}/enrichment_summary_top5_per_analysis.csv", index=False)
        print(f"\nSummary table written")

    print("\n=== g:Profiler enrichment complete ===")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--lasso_dir",  required=True)
    p.add_argument("--de_csv",     required=True)
    p.add_argument("--grn_csv",    required=True)
    p.add_argument("--outdir",     required=True)
    p.add_argument("--tfs",        nargs="+", default=["BACH2", "ATF4", "XBP1"])
    p.add_argument("--n_de",       type=int, default=200)
    p.add_argument("--top_terms",  type=int, default=15)
    args = p.parse_args()
    main(args.lasso_dir, args.de_csv, args.grn_csv,
         args.outdir, args.tfs, args.n_de, args.top_terms)