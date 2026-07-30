#!/usr/bin/env python
"""
04_aggregate_grn_scores.py

Pools the GRNBoost2 'adjacencies.tsv' (TF, target, importance) from every
donor x cell_type pySCENIC run, finds the most recurrent TF-target pairs
(i.e. pairs that show up as edges across the most runs), and visualises
their importance-score distribution per cell type with violin+inset-boxplot.

Usage:
    python 04_aggregate_grn_scores.py --pyscenic_dir ./pyscenic_out \
                                        --outdir ./grn_analysis_out \
                                        --top_n 20
"""
import argparse
import glob
import os
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def load_all_adjacencies(pyscenic_dir):
    records = []
    for adj_path in glob.glob(f"{pyscenic_dir}/*/adjacencies.tsv"):
        run_name = os.path.basename(os.path.dirname(adj_path))   # e.g. BM1_Memory_B
        m = re.match(r"(BM\d+)_(Memory_B|Plasma)", run_name)
        if not m:
            continue
        donor, cell_type = m.group(1), m.group(2)
        df = pd.read_csv(adj_path, sep="\t")  # columns: TF, target, importance
        df["donor"] = donor
        df["cell_type"] = cell_type
        df["pair"] = df["TF"] + "->" + df["target"]
        records.append(df)
    return pd.concat(records, ignore_index=True)


def main(pyscenic_dir, outdir, top_n):
    os.makedirs(outdir, exist_ok=True)
    all_adj = load_all_adjacencies(pyscenic_dir)
    all_adj.to_parquet(f"{outdir}/all_adjacencies.parquet")

    # how many distinct runs (donor x cell_type) exist a given pair in
    n_runs_per_pair = all_adj.groupby(["pair", "cell_type"])["donor"].nunique().reset_index(name="n_donors")
    most_common = (n_runs_per_pair.sort_values(["cell_type", "n_donors"], ascending=[True, False])
                   .groupby("cell_type").head(top_n))
    most_common.to_csv(f"{outdir}/most_common_pairs.csv", index=False)

    # summary stats (mean, sd) of importance for the most common pairs
    plot_df = all_adj.merge(most_common[["pair", "cell_type"]], on=["pair", "cell_type"])
    summary = plot_df.groupby(["pair", "cell_type"])["importance"].agg(["mean", "std", "count"]).reset_index()
    summary.to_csv(f"{outdir}/most_common_pairs_importance_summary.csv", index=False)

    # --- violin + inset boxplot, one panel per cell type ---
    for cell_type in plot_df.cell_type.unique():
        sub = plot_df[plot_df.cell_type == cell_type]
        order = (sub.groupby("pair")["importance"].mean()
                 .sort_values(ascending=False).index)

        fig, ax = plt.subplots(figsize=(max(8, 0.5 * len(order)), 6))
        sns.violinplot(data=sub, x="pair", y="importance", order=order,
                        inner=None, cut=0, ax=ax, color="lightsteelblue")
        sns.boxplot(data=sub, x="pair", y="importance", order=order,
                     width=0.15, showcaps=True, boxprops={"zorder": 2, "facecolor": "white"},
                     whiskerprops={"zorder": 2}, ax=ax, fliersize=1)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_title(f"Top {top_n} most recurrent TF-target pairs - {cell_type}")
        ax.set_ylabel("GRNBoost2 importance score")
        fig.tight_layout()
        fig.savefig(f"{outdir}/violin_box_{cell_type}.png", dpi=200)
        plt.close(fig)

    print("Done. See most_common_pairs.csv and violin_box_*.png")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--pyscenic_dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--top_n", type=int, default=20)
    args = p.parse_args()
    main(args.pyscenic_dir, args.outdir, args.top_n)
