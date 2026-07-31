#!/usr/bin/env python
"""
04_aggregate_grn_scores.py

This script pools the GRNBoost2 adjacencies.tsv's from every
donor x cell_type pySCENIC run, finds the most recurrent TF-target pairs and visualises
their importance-score distribution per cell type with violin+inset-boxplot.

Usage:
    python 04_aggregate_grn_scores.py --pyscenic_dir /path/to/pyscenic_out \
                                       --outdir /path/to/grn_analysis_out \
                                       --top_n 20
"""
import argparse
import glob
import os
import re
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")  
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    'font.size': 18,
    'axes.titlesize': 20,
    'axes.labelsize': 18,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
})


def load_all_adjacencies(pyscenic_dir):
    records = []
    for adj_path in sorted(glob.glob(f"{pyscenic_dir}/*/adjacencies.tsv")):
        run_name = os.path.basename(os.path.dirname(adj_path))
        m = re.match(r"(BM\d+)_(Memory_B|Plasma)", run_name)
        if not m:
            continue
        donor, cell_type = m.group(1), m.group(2)
        df = pd.read_csv(adj_path, sep="\t")
        df.columns = df.columns.str.strip()
        df["donor"] = donor
        df["cell_type"] = cell_type
        df["pair"] = df["TF"].astype(str) + "->" + df["target"].astype(str)
        records.append(df)
        print(f"Loaded {run_name}: {len(df)} edges")
    if not records:
        raise ValueError(f"No adjacencies.tsv files found under {pyscenic_dir}")
    return pd.concat(records, ignore_index=True)


def main(pyscenic_dir, outdir, top_n):
    os.makedirs(outdir, exist_ok=True)

    print("Loading all adjacencies...")
    all_adj = load_all_adjacencies(pyscenic_dir)
    print(f"Total edges loaded: {len(all_adj)}")

    # save pooled table
    all_adj.to_csv(f"{outdir}/all_adjacencies.csv", index=False)

    # rank pairs by how many donors they appear in, per cell type
    n_donors_per_pair = all_adj.groupby(["pair", "cell_type"])["donor"].nunique().reset_index(name="n_donors")
    most_common = (n_donors_per_pair.sort_values(["cell_type", "n_donors"], ascending=[True, False]).groupby("cell_type").head(top_n))
    most_common.to_csv(f"{outdir}/most_common_pairs.csv", index=False)
    print(f"Most common pairs saved.")

    # mean/sd importance for top pairs
    plot_df = all_adj.merge(most_common[["pair", "cell_type"]], on=["pair", "cell_type"])
    summary = plot_df.groupby(["pair", "cell_type"])["importance"].agg(["mean", "std", "count"]).reset_index()
    summary.to_csv(f"{outdir}/most_common_pairs_importance_summary.csv", index=False)
    print(f"Summary stats saved.")

    # violin + inset boxplot per cell type
    for cell_type in sorted(plot_df.cell_type.unique()):
        sub = plot_df[plot_df.cell_type == cell_type].copy()
        order = list(sub.groupby("pair")["importance"].mean().sort_values(ascending=False).index)

        fig, ax = plt.subplots(figsize=(max(10, 0.8 * len(order)), 7))

        sns.violinplot(data=sub, x="pair", y="importance", order=order, inner=None, cut=0, ax=ax, color="lightsteelblue", alpha=0.8)
        sns.boxplot(data=sub, x="pair", y="importance", order=order, width=0.12, showcaps=True, ax=ax, fliersize=1, boxprops={"zorder": 2, "facecolor": "white"}, whiskerprops={"zorder": 2}, medianprops={"zorder": 3, "color": "red"})

        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=16)
        ax.set_title(f"Top {top_n} most recurrent TF-target pairs — {cell_type}", fontsize=18)
        ax.set_ylabel("GRNBoost2 importance score", fontsize=16)
        ax.set_xlabel("TF -> target pair", fontsize=16)
        fig.tight_layout()
        out_path = f"{outdir}/violin_box_{cell_type}.png"
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"Plot saved: {out_path}")

    print("Done.")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--pyscenic_dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--top_n", type=int, default=20)
    args = p.parse_args()
    main(args.pyscenic_dir, args.outdir, args.top_n)
