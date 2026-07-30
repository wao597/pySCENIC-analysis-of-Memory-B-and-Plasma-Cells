#!/usr/bin/env python
"""
08_tf_network_focused.py

This script builds two focused TF regulon networks: BACH2 in Memory B cells and ATF4 in Plasma cells by using the most recurrent TF-target pairs from most_common_pairs.csv, with edge width scaled by mean importance score across donors from the pooled adjacencies.


Output:
  - BACH2_Memory_B_network.png
  - ATF4_Plasma_network.png
  - BACH2_ATF4_combined_panel.png  

Usage:
    python 08_tf_network_focused.py \
        --pairs_csv /path/to/grn_analysis_out/most_common_pairs.csv \
        --adj_dir   /path/to/pyscenic_out \
        --outdir    /path/to/network_out \
        --top_targets 10 \
        --min_donors 6
"""
import argparse
import os
import glob
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import networkx as nx


# visual colour specifics  
BACH2_COLOUR   = "#1565C0"   # blue   — Memory B
ATF4_COLOUR    = "#C62828"   # red    — Plasa
TF_NODE_SIZE   = 4000
TARGET_SIZE    = 900
TARGET_COLOUR  = "#F5F5F5"
TARGET_OUTLINE = "#90A4AE"
BACKGROUND     = "white"


# data loading 
def load_mean_importance(adj_dir, cell_type, tf):
    """
    Pool adjacencies for one cell_type across all 8 donors.
    Returns mean importance score per TF-target pair.
    """
    records = []
    for path in glob.glob(f"{adj_dir}/*_{cell_type}/adjacencies.tsv"):
        df = pd.read_csv(path, sep="\t")
        df.columns = df.columns.str.strip()
        df = df[df["TF"] == tf]
        if len(df) > 0:
            records.append(df)
    if not records:
        print(f"  WARNING: no adjacencies found for {tf} in {cell_type}",
              flush=True)
        return pd.DataFrame(columns=["target", "mean_importance"])
    pooled = pd.concat(records, ignore_index=True)
    summary = (pooled.groupby("target")["importance"]
               .agg(mean_importance="mean", n_donors="count")
               .reset_index())
    return summary


def load_recurrent_targets(pairs_csv, adj_dir, tf, cell_type, top_n, min_donors):

    pairs = pd.read_csv(pairs_csv)
    if "pair" in pairs.columns:
        pairs[["TF", "target"]] = pairs["pair"].str.split("->", expand=True)
    sub = pairs[(pairs["TF"] == tf) &
                (pairs["cell_type"] == cell_type) &
                (pairs["n_donors"] >= min_donors)].copy()

    # Deal with ties by using mean importance score across donors
    importance_df = load_mean_importance(adj_dir, cell_type, tf)
    if len(importance_df) > 0:
        sub = sub.merge(importance_df[["target", "mean_importance"]],
                        on="target", how="left")
        sub = sub.sort_values(["n_donors", "mean_importance"],
                               ascending=[False, False])
    else:
        sub = sub.sort_values("n_donors", ascending=False)

    sub = sub.head(top_n)
    print(f"  {tf} → {cell_type}: {len(sub)} targets "
          f"(≥{min_donors}/8 donors, top {top_n} by importance)", flush=True)
    return sub


# Network construction 
def build_regulon_graph(targets_df, importance_df, tf):
    """Build a star-topology DiGraph: TF → targets."""
    G = nx.DiGraph()
    G.add_node(tf, node_type="TF")

    for _, row in targets_df.iterrows():
        target = row["target"]
        n_donors = row["n_donors"]

        imp_row = importance_df[importance_df["target"] == target]
        mean_imp = float(imp_row["mean_importance"].iloc[0]) \
            if len(imp_row) > 0 else 1.0

        G.add_node(target, node_type="target")
        G.add_edge(tf, target,
                   weight=mean_imp,
                   n_donors=n_donors)
    return G


# Drawing 
def draw_regulon(G, tf, tf_colour, cell_type, ax, title, top_n):
    """Draw a single TF regulon network onto ax."""
    if G.number_of_nodes() == 0:
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, fontsize=12)
        ax.set_title(title)
        ax.axis("off")
        return

    target_nodes = [n for n in G.nodes() if n != tf]
    n_targets = len(target_nodes)

    # circular layout - TF at centre, targets around it
    pos = {tf: np.array([0.0, 0.0])}
    angles = np.linspace(0, 2 * np.pi, n_targets, endpoint=False)
    radius = 1.6
    for node, angle in zip(target_nodes, angles):
        pos[node] = np.array([radius * np.cos(angle),
                               radius * np.sin(angle)])

    # edge widths scaled by importance
    edges = list(G.edges(data=True))
    weights = [d["weight"] for _, _, d in edges]
    max_w = max(weights) if weights else 1.0
    min_w = min(weights) if weights else 0.0
    rng = max_w - min_w if max_w > min_w else 1.0
    edge_widths = [1.0 + 5.0 * (w - min_w) / rng for w in weights]
    edge_alphas = [0.4 + 0.5 * (w - min_w) / rng for w in weights]

    # colour edges by importance
    norm_weights = [(w - min_w) / rng for w in weights]
    base = matplotlib.colors.to_rgb(tf_colour)
    edge_colours = [tuple(c * nw + 0.9 * (1 - nw) for c in base)
                    for nw in norm_weights]

    # draw target nodes
    nx.draw_networkx_nodes(G, pos, ax=ax,
                           nodelist=target_nodes,
                           node_color=TARGET_COLOUR,
                           node_size=TARGET_SIZE,
                           edgecolors=TARGET_OUTLINE,
                           linewidths=1.2)

    # draw TF node on top
    nx.draw_networkx_nodes(G, pos, ax=ax,
                           nodelist=[tf],
                           node_color=tf_colour,
                           node_size=TF_NODE_SIZE,
                           edgecolors="white",
                           linewidths=2.5)

    # draw edges
    for (u, v, d), width, colour, alpha in zip(edges, edge_widths,
                                                 edge_colours, edge_alphas):
        nx.draw_networkx_edges(
            G, pos, ax=ax,
            edgelist=[(u, v)],
            width=width,
            edge_color=[colour],
            alpha=alpha,
            arrows=True,
            arrowsize=10,
            arrowstyle="-|>",
            min_source_margin=40,
            min_target_margin=18
        )

    # TF label 
    nx.draw_networkx_labels(G, pos, ax=ax,
                            labels={tf: tf},
                            font_size=16,
                            font_weight="bold",
                            font_color="white")

    # target labels 
    nx.draw_networkx_labels(G, pos, ax=ax,
                            labels={n: n for n in target_nodes},
                            font_size=10,
                            font_color="#263238",
                            font_weight="bold")

    # colourbar for edge importance
    sm = cm.ScalarMappable(
        cmap=matplotlib.colors.LinearSegmentedColormap.from_list(
            "imp", ["#ECEFF1", tf_colour]),
        norm=plt.Normalize(vmin=min_w, vmax=max_w)
    )
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.35, pad=0.02, aspect=15)
    cbar.set_label("Mean GRNBoost2\nimportance score", fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    # legend
    tf_patch = mpatches.Patch(color=tf_colour, label=f"{tf} (TF)")
    tg_patch = mpatches.Patch(facecolor=TARGET_COLOUR,
                               edgecolor=TARGET_OUTLINE,
                               label="Target gene")
    ax.legend(handles=[tf_patch, tg_patch],
              loc="upper right", fontsize=11, framealpha=0.85)

    ax.set_title(title, fontsize=13, fontweight="bold", pad=14)
    ax.axis("off")


# Main Network 
def make_one_network(pairs_csv, adj_dir, outdir,
                     tf, cell_type, tf_colour,
                     top_targets, min_donors):
    print(f"\n=== {tf} — {cell_type} ===", flush=True)

    targets_df    = load_recurrent_targets(
        pairs_csv, adj_dir, tf, cell_type, top_targets, min_donors)
    importance_df = load_mean_importance(adj_dir, cell_type, tf)
    G             = build_regulon_graph(targets_df, importance_df, tf)

    title = (f"{tf} regulon — {cell_type.replace('_', ' ')} cells\n"
             f"Top {len(targets_df)} targets recurrent in "
             f"≥{min_donors}/8 donors")

    fig, ax = plt.subplots(figsize=(14, 13), facecolor=BACKGROUND)
    ax.set_facecolor(BACKGROUND)
    draw_regulon(G, tf, tf_colour, cell_type, ax, title, top_targets)
    fig.tight_layout()
    out = f"{outdir}/{tf}_{cell_type}_network.png"
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor=BACKGROUND)
    plt.close(fig)
    print(f"  Saved: {out}", flush=True)

    return G


def main(pairs_csv, adj_dir, outdir, top_targets, min_donors):
    os.makedirs(outdir, exist_ok=True)

    # individual networks
    make_one_network(
        pairs_csv, adj_dir, outdir,
        tf="BACH2", cell_type="Memory_B",
        tf_colour=BACH2_COLOUR,
        top_targets=top_targets, min_donors=min_donors)

    make_one_network(
        pairs_csv, adj_dir, outdir,
        tf="ATF4", cell_type="Plasma",
        tf_colour=ATF4_COLOUR,
        top_targets=top_targets, min_donors=min_donors)

    # combined visualisation of networks
    print("\n=== Combined panel ===", flush=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(26, 13),
                                    facecolor=BACKGROUND)
    for ax in [ax1, ax2]:
        ax.set_facecolor(BACKGROUND)

    t_bach2 = load_recurrent_targets(
        pairs_csv, adj_dir, "BACH2", "Memory_B", top_targets, min_donors)
    i_bach2 = load_mean_importance(adj_dir, "Memory_B", "BACH2")
    G_b = build_regulon_graph(t_bach2, i_bach2, "BACH2")
    draw_regulon(G_b, "BACH2", BACH2_COLOUR, "Memory_B", ax1,
                 f"BACH2 regulon — Memory B cells\n"
                 f"(top {len(t_bach2)} targets, ≥{min_donors}/8 donors)",
                 top_targets)

    t_atf4 = load_recurrent_targets(
        pairs_csv, adj_dir, "ATF4", "Plasma", top_targets, min_donors)
    i_atf4 = load_mean_importance(adj_dir, "Plasma", "ATF4")
    G_a = build_regulon_graph(t_atf4, i_atf4, "ATF4")
    draw_regulon(G_a, "ATF4", ATF4_COLOUR, "Plasma", ax2,
                 f"ATF4 regulon — Plasma cells\n"
                 f"(top {len(t_atf4)} targets, ≥{min_donors}/8 donors)",
                 top_targets)

    fig.suptitle(
        "Recurrent TF Regulon Networks: Memory B cells vs Plasma cells\n"
        "Edge width and colour intensity = mean GRNBoost2 importance score",
        fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    out = f"{outdir}/BACH2_ATF4_combined_panel.png"
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor=BACKGROUND)
    plt.close(fig)
    print(f"  Saved: {out}", flush=True)
    print("\n=== Done ===")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--pairs_csv",   required=True)
    p.add_argument("--adj_dir",     required=True)
    p.add_argument("--outdir",      required=True)
    p.add_argument("--top_targets", type=int, default=10,
                   help="Max target genes to show per TF (default 10)")
    p.add_argument("--min_donors",  type=int, default=6,
                   help="Min donors edge must appear in (default 6)")
    args = p.parse_args()
    main(args.pairs_csv, args.adj_dir, args.outdir,
         args.top_targets, args.min_donors)
