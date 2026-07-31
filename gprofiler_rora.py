#!/usr/bin/env python
"""
gprofiler_rora.py
This script runs grprofiler on the lasso derived RORA terms
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
    cleaned = []
    for g in gene_list:
        if pd.isna(g) or g is None:
            continue
        g_str = str(g).strip()
        if g_str == "" or g_str.lower() == "nan" or g_str == "None":
            continue
        if g_str[0].isdigit():
            continue
        cleaned.append(g_str)
    return cleaned

def convert_ensg_to_symbol(genes):
    if not genes:
        return []
    
    ensg_list = [g for g in genes if str(g).startswith('ENSG')]
    symbol_list = [g for g in genes if not str(g).startswith('ENSG')]
    
    if not ensg_list:
        return genes
    
    try:
        gp = GProfiler(return_dataframe=True)
        converted = gp.convert(organism="hsapiens", query=ensg_list)
        
        if converted is not None and not converted.empty:
            mapping = dict(zip(converted['incoming'], converted['converted']))
            converted_symbols = [mapping.get(g, g) for g in ensg_list]
            all_genes = symbol_list + converted_symbols
            all_genes = [g for g in all_genes if not str(g).startswith('ENSG')]
            return all_genes
        else:
            return genes
    except Exception as e:
        print(f"Conversion error: {e}")
        return genes

def run_gprofiler(gene_list, label, outdir, top_terms=20):
    gene_list = clean_gene_list(gene_list)
    if len(gene_list) == 0:
        print(f"Empty gene list for {label}")
        return None
    
    gene_list = convert_ensg_to_symbol(gene_list)
    gene_list = list(dict.fromkeys(gene_list))
    
    if len(gene_list) == 0:
        print(f"No valid genes after conversion for {label}")
        return None
    
    print(f"Running g:Profiler: {label} ({len(gene_list)} genes)")
    
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
        print(f"Error running g:Profiler: {e}")
        return None
    
    safe_label = label.replace(" ", "_").replace("/", "_").replace("-", "_")

    if results is None or results.empty:
        print(f"No significant terms found for {label}")
        return None
    
    results.to_csv(f"{outdir}/{safe_label}_gprofiler_results.csv", index=False)
    print(f"Found {len(results)} significant terms")

    
    display_label = label.replace("_", " ")
    
    # Visualisation
    plot_df = results.sort_values("p_value").head(top_terms).copy()
    plot_df["gene_ratio"] = plot_df["intersection_size"] / plot_df["query_size"]
    plot_df["-log10_padj"] = -np.log10(plot_df["p_value"].clip(lower=1e-300))
    plot_df["name"] = plot_df["name"].str[:60]

    # sort by significance so most significant comes first 
    plot_df = plot_df.sort_values("-log10_padj", ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(10, max(6, 0.5 * len(plot_df))),
                           facecolor="white")

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
    ax.invert_yaxis()  #
    ax.set_xlabel("Gene Ratio (Intersection / Query Size)", fontsize=11, fontweight="bold")
    ax.set_title(f"Top Enriched Terms: {display_label}", fontsize=12, fontweight="bold", pad=10)

    cbar = fig.colorbar(scatter, ax=ax, shrink=0.4, pad=0.02, aspect=18)
    cbar.set_label("-log10(adjusted p-value)\nDarker = more significant", fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    for size in [5, 20, 50]:
        ax.scatter([], [], s=size * 15, c="grey", alpha=0.6, label=f"n={size} genes")
    ax.legend(title="Intersection size", title_fontsize=9,
              loc="lower right", fontsize=8, framealpha=0.85)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#ECEFF1", linewidth=0.7)

    fig.tight_layout()
    fig.savefig(f"{outdir}/{safe_label}_dotplot.png", dpi=300, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    
    return results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lasso_dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--cell_type", default="Memory_B", choices=["Memory_B", "Plasma"])
    parser.add_argument("--top_terms", type=int, default=20)
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    tf = "RORA"
    coef_path = f"{args.lasso_dir}/{tf}/{args.cell_type}_{tf}_lasso_coefs.csv"
    
    if not os.path.exists(coef_path):
        print(f"Error: {coef_path} not found")
        return
    
    coefs = pd.read_csv(coef_path, index_col=0, header=None)
    coefs.columns = ["coefficient"]
    coefs["coefficient"] = pd.to_numeric(coefs["coefficient"], errors="coerce")
    coefs = coefs.dropna(subset=["coefficient"])
    coefs = coefs[coefs["coefficient"] != 0]
    
    genes = coefs.index.tolist()
    label = f"RORA{args.cell_type}Lasso"
    
    results = run_gprofiler(genes, label, args.outdir, args.top_terms)
    
    if results is not None:
        print(f"\nTop 10 terms for RORA {args.cell_type}:")
        print(results[["name", "source", "p_value", "intersection_size"]].head(10).to_string(index=False))
    
    print(f"\nResults saved to: {args.outdir}")

if __name__ == "__main__":
    main()
