#!/usr/bin/env python3
"""
run_18s_pipeline.py
===================
A cohesive end-to-end pipeline for 18S rRNA extraction, classification, 
diversity visualization, and target contig extraction.
"""

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

# ---------------------------------------------------------------------------
# Global Plotting Configuration
# ---------------------------------------------------------------------------
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "axes.edgecolor": "#333333",
    "axes.labelcolor": "#222222",
    "text.color": "#222222",
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

ALL_RANKS = ["Domain", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"]
BASE_PALETTE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2",
    "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
    "#1B7837", "#B2182B", "#2166AC", "#D6604D", "#5AAE61",
]
TAXON_RE = re.compile(r"^\s*(.*?)\s*\(([-+]?[0-9]*\.?[0-9]+)\)\s*$")


# ---------------------------------------------------------------------------
# Pipeline Execution Functions
# ---------------------------------------------------------------------------
def run_command(cmd, stdout_path=None):
    """Helper to run shell commands via subprocess."""
    try:
        if stdout_path:
            with open(stdout_path, 'w') as f:
                subprocess.run(cmd, stdout=f, check=True, text=True)
        else:
            subprocess.run(cmd, check=True, text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error executing command: {' '.join(cmd)}\n{e.stderr}")
    except FileNotFoundError:
        sys.exit(f"Error: Required tool '{cmd[0]}' is not installed or not in PATH.")

def filter_18s_gff(in_gff: Path, out_gff: Path) -> int:
    """Replaces awk to filter barrnap GFF output for 18S only."""
    count = 0
    with open(in_gff, 'r') as infile, open(out_gff, 'w') as outfile:
        for line in infile:
            if line.startswith('#'): continue
            parts = line.split('\t')
            if len(parts) > 8 and parts[2] == "rRNA" and "Name=18S_rRNA" in parts[8]:
                outfile.write(line)
                count += 1
    return count

def parse_sintax_to_csv(input_file: Path, output_file: Path):
    """Parses raw VSEARCH SINTAX tabbed output into a clean taxonomy CSV."""
    if not input_file.exists() or input_file.stat().st_size == 0:
        print("Notice: No classification results generated.")
        return False

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["Contig_ID", "Domain", "Supergroup", "Division", "Class", 
                         "Order", "Family", "Genus", "Species", "Full_Prediction_String"])
        
        for line in infile:
            if not line.strip(): continue
            parts = line.strip().split('\t')
            contig_id = parts[0]
            prediction_string = parts[1] if len(parts) > 1 else ""
            
            tax_dict = {k: 'Unclassified' for k in ['d', 'k', 'p', 'c', 'o', 'f', 'g', 's']}
            
            if prediction_string:
                for taxon in prediction_string.split(','):
                    if ':' in taxon:
                        rank, name = taxon.split(':', 1)
                        if rank in tax_dict:
                            tax_dict[rank] = name
            
            writer.writerow([
                contig_id, tax_dict['d'], tax_dict['k'], tax_dict['p'],
                tax_dict['c'], tax_dict['o'], tax_dict['f'],
                tax_dict['g'], tax_dict['s'], prediction_string
            ])
    return True


# ---------------------------------------------------------------------------
# Plotting & Diversity Functions
# ---------------------------------------------------------------------------
def parse_taxon_field(value):
    if pd.isna(value): return pd.NA, np.nan
    s = str(value).strip()
    m = TAXON_RE.match(s)
    if m: return m.group(1), float(m.group(2))
    return s, np.nan

def load_and_clean(csv_path: Path, ranks=ALL_RANKS, drop_bacteria=True, drop_mito=True):
    df = pd.read_csv(csv_path)
    present_ranks = [r for r in ranks if r in df.columns]
    
    for rank in present_ranks:
        parsed = df[rank].apply(parse_taxon_field)
        df[f"{rank}_name"] = parsed.apply(lambda t: t[0])
        df[f"{rank}_confidence"] = parsed.apply(lambda t: t[1])

    if drop_bacteria and "Supergroup_name" in df.columns:
        is_bacteria = df["Supergroup_name"].astype(str).str.strip().str.lower() == "bacteria"
        df = df.loc[~is_bacteria].copy()

    if drop_mito and "Supergroup" in df.columns:
        mito_mask = df["Supergroup"].astype(str).str.contains("mito", case=False, na=False)
        df = df.loc[~mito_mask].copy()

    return df.reset_index(drop=True)

def diversity_stats(counts: pd.Series) -> dict:
    counts = counts[counts > 0]
    n = counts.sum()
    p = counts / n
    richness = len(counts)
    return {
        "n_contigs": int(n),
        "richness": richness,
        "shannon_H": -np.sum(p * np.log(p)) if richness > 0 else 0.0,
        "simpson_1_minus_D": 1 - np.sum(p ** 2) if richness > 0 else 0.0,
        "pielou_evenness": (-np.sum(p * np.log(p))) / np.log(richness) if richness > 1 else np.nan,
    }

def _palette(n):
    if n <= len(BASE_PALETTE): return BASE_PALETTE[:n]
    return (BASE_PALETTE * int(np.ceil(n / len(BASE_PALETTE))))[:n]

def _counts_for_rank(df: pd.DataFrame, rank: str, top_n: int) -> pd.Series:
    counts = df[f"{rank}_name"].dropna().value_counts()
    if len(counts) > top_n:
        top = counts.iloc[:top_n]
        counts = pd.concat([top, pd.Series({"Other": counts.iloc[top_n:].sum()})])
    return counts

def plot_rank_barh(df, rank, top_n, out_dir: Path, dpi):
    counts = _counts_for_rank(df, rank, top_n)
    if counts.empty: return None

    colors = _palette(len(counts))
    if "Other" in counts.index: colors[list(counts.index).index("Other")] = "#B0B0B0"

    fig, ax = plt.subplots(figsize=(8, max(2.5, 0.42 * len(counts) + 1.2)))
    order = counts.sort_values(ascending=True)
    
    bars = ax.barh(order.index, order.values, color=[colors[list(counts.index).index(n)] for n in order.index],
                   edgecolor="white", linewidth=0.6, height=0.72)

    total = order.values.sum()
    for bar, val in zip(bars, order.values):
        ax.text(bar.get_width() + total * 0.01, bar.get_y() + bar.get_height() / 2,
                f"{val:,} ({100*val/total:.1f}%)", va="center", ha="left", fontsize=8.5, color="#333333")

    ax.set_xlabel("Number of contigs")
    ax.set_title(f"18S taxonomic diversity — {rank}\n(n={total:,} contigs)", fontsize=12, fontweight="bold", loc="left")
    ax.set_xlim(0, order.values.max() * 1.22)
    ax.grid(axis="x", color="#E0E0E0", linewidth=0.7)
    ax.set_axisbelow(True)
    fig.tight_layout()

    stem = out_dir / f"barh_{rank.lower()}"
    fig.savefig(f"{stem}.pdf", dpi=dpi, bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

def plot_rank_donut(df, rank, top_n, out_dir: Path, dpi):
    counts = _counts_for_rank(df, rank, top_n)
    if counts.empty: return None

    counts = counts.sort_values(ascending=False)
    colors = _palette(len(counts))
    if "Other" in counts.index: colors[list(counts.index).index("Other")] = "#B0B0B0"

    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    total = counts.sum()
    wedges, _ = ax.pie(counts.values, colors=colors, startangle=90, counterclock=False,
                       wedgeprops=dict(width=0.42, edgecolor="white", linewidth=1.2))

    for wedge, name, val in zip(wedges, counts.index, counts.values):
        frac = val / total
        if frac >= 0.025:
            ang = np.deg2rad((wedge.theta1 + wedge.theta2) / 2)
            ax.annotate(f"{name}\n{100*frac:.1f}%", xy=(0.79 * np.cos(ang), 0.79 * np.sin(ang)),
                        xytext=(1.18 * np.cos(ang), 1.05 * np.sin(ang)), ha="center", va="center", fontsize=8.5,
                        arrowprops=dict(arrowstyle="-", color="#999999", lw=0.8))

    ax.text(0, 0, f"{total:,}\ncontigs", ha="center", va="center", fontsize=13, fontweight="bold")
    ax.set_title(f"Composition by {rank}", fontsize=12, fontweight="bold")
    ax.set_aspect("equal")
    fig.tight_layout()

    stem = out_dir / f"donut_{rank.lower()}"
    fig.savefig(f"{stem}.pdf", dpi=dpi, bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

def plot_overview(df, ranks, top_n, out_dir: Path, dpi, stats_table: pd.DataFrame):
    ncols = 3
    nrows = int(np.ceil((len(ranks) + 1) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.2 * ncols, 4.4 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for ax, rank in zip(axes, ranks):
        counts = _counts_for_rank(df, rank, top_n=8)
        if counts.empty:
            ax.axis("off")
            continue
        counts = counts.sort_values(ascending=True)
        colors = _palette(len(counts))
        if "Other" in counts.index: colors[list(counts.index).index("Other")] = "#B0B0B0"
        
        ax.barh(counts.index, counts.values, color=colors, edgecolor="white", linewidth=0.5)
        ax.set_title(rank, fontsize=11, fontweight="bold", loc="left")
        ax.tick_params(labelsize=8)
        ax.grid(axis="x", color="#E5E5E5", linewidth=0.6)
        ax.set_axisbelow(True)

    stats_ax = axes[len(ranks)]
    stats_ax.axis("off")
    cell_text = [[row.Index, f"{row.richness}", f"{row.shannon_H:.2f}", 
                  f"{row.simpson_1_minus_D:.2f}", f"{row.pielou_evenness:.2f}" if not np.isnan(row.pielou_evenness) else "–"] 
                 for row in stats_table.itertuples()]
    
    table = stats_ax.table(cellText=cell_text, colLabels=["Rank", "Richness", "Shannon H'", "Simpson", "Pielou"], 
                           loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.6)
    stats_ax.set_title("Diversity indices", fontsize=11, fontweight="bold", loc="left")

    for ax in axes[len(ranks) + 1:]: ax.axis("off")

    fig.suptitle("18S Taxonomic Diversity — Overview", fontsize=15, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(out_dir / "overview_all_ranks.pdf", dpi=dpi, bbox_inches="tight")
    fig.savefig(out_dir / "overview_all_ranks.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

# ---------------------------------------------------------------------------
# Contig Extraction Function
# ---------------------------------------------------------------------------
def extract_full_contigs(assembly_path: Path, filtered_df: pd.DataFrame, out_fasta: Path):
    """
    Parses the final filtered DataFrame to map 18S sequences back to their 
    parent contigs in the original assembly, and saves them to a new FASTA.
    """
    clean_ids = set()
    for cid in filtered_df["Contig_ID"]:
        # bedtools getfasta appends coordinates like ':start-end' to the ID. 
        # We strip this using regex to get the raw parent contig ID.
        base_id = re.sub(r':\d+-\d+.*$', '', str(cid))
        clean_ids.add(base_id)
        
    print(f"      Identified {len(clean_ids)} unique eukaryotic contigs for extraction.")
    
    extracted_count = 0
    with open(out_fasta, "w") as out_f:
        # Bio.SeqIO handles the raw assembly parsing
        for record in SeqIO.parse(assembly_path, "fasta"):
            if record.id in clean_ids:
                SeqIO.write(record, out_f, "fasta")
                extracted_count += 1
                
    print(f"      Successfully saved {extracted_count} full contigs to {out_fasta.name}")


# ---------------------------------------------------------------------------
# Main Orchestrator
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="End-to-end 18S extraction, classification, and plotting pipeline.")
    parser.add_argument("-a", "--assembly", type=Path, required=True, help="Path to input long-read assembly FASTA")
    parser.add_argument("-d", "--db", type=Path, required=True, help="Path to UTAX formatted PR2 database FASTA")
    parser.add_argument("-o", "--out-dir", type=Path, default=Path("18s_output"), help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of CPU threads")
    parser.add_argument("-c", "--cutoff", type=float, default=0.8, help="SINTAX confidence cutoff")
    
    # Plotting arguments
    parser.add_argument("--ranks", nargs="+", default=["Division", "Class", "Order", "Family", "Genus"], choices=ALL_RANKS)
    parser.add_argument("--top-n", type=int, default=15, help="Max taxa to show individually per rank")
    parser.add_argument("--keep-bacteria", action="store_true", help="Keep Bacterial hits")
    parser.add_argument("--keep-mito", action="store_true", help="Keep Mitochondrial hits")
    parser.add_argument("--dpi", type=int, default=300, help="DPI for PNG plots")
    parser.add_argument("--no-donuts", action="store_true", help="Skip donut charts")
    parser.add_argument("--no-overview", action="store_true", help="Skip overview figure")
    
    args = parser.parse_args()

    if not args.assembly.exists(): sys.exit(f"Error: Assembly {args.assembly} not found.")
    if not args.db.exists(): sys.exit(f"Error: Database {args.db} not found.")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    
    # File Paths
    gff_out = args.out_dir / "rRNA_predictions.gff"
    gff_18s = args.out_dir / "18S_only.gff"
    fasta_18s = args.out_dir / "18S_sequences.fasta"
    sintax_out = args.out_dir / "18S_sintax_output.txt"
    csv_out = args.out_dir / "18S_taxonomy_results.csv"

    # Step 1: Barrnap
    print("[1/6] Running Barrnap to predict eukaryotic rRNA...")
    run_command(["barrnap", "--kingdom", "euk", "--threads", str(args.threads), str(args.assembly)], stdout_path=gff_out)

    # Step 2: Filter GFF
    print("[2/6] Filtering predictions for 18S rRNA sequences...")
    num_18s = filter_18s_gff(gff_out, gff_18s)
    if num_18s == 0:
        print("Notice: No 18S sequences found. Exiting gracefully.")
        sys.exit(0)

    # Step 3: Bedtools Extraction
    print(f"[3/6] Extracting {num_18s} 18S sequences...")
    run_command(["bedtools", "getfasta", "-fi", str(args.assembly), "-bed", str(gff_18s), "-fo", str(fasta_18s)])

    # Step 4: VSEARCH Classification
    print("[4/6] Classifying sequences against PR2 database via VSEARCH...")
    run_command(["vsearch", "--sintax", str(fasta_18s), "--db", str(args.db), 
                 "--tabbedout", str(sintax_out), "--sintax_cutoff", str(args.cutoff), "--threads", str(args.threads)])

    # Step 5: Convert and Plot
    print("[5/6] Parsing output and generating diversity plots...")
    if parse_sintax_to_csv(sintax_out, csv_out):
        df = load_and_clean(csv_out, ranks=ALL_RANKS, drop_bacteria=not args.keep_bacteria, drop_mito=not args.keep_mito)
        if df.empty:
            print("Notice: No valid eukaryotic contigs remain after filtering.")
            sys.exit(0)
            
        df.to_csv(args.out_dir / "filtered_taxonomy_table.csv", index=False)
        
        stats_rows = {}
        for rank in args.ranks:
            col = f"{rank}_name"
            if col not in df.columns: continue
            
            stats_rows[rank] = diversity_stats(df[col].dropna().value_counts())
            plot_rank_barh(df, rank, args.top_n, args.out_dir, args.dpi)
            if not args.no_donuts: plot_rank_donut(df, rank, args.top_n, args.out_dir, args.dpi)
            
        stats_table = pd.DataFrame.from_dict(stats_rows, orient="index")
        stats_table.index.name = "Rank"
        stats_table.to_csv(args.out_dir / "diversity_statistics.csv")
        
        if not args.no_overview and stats_rows:
            plot_overview(df, list(stats_rows.keys()), args.top_n, args.out_dir, args.dpi, stats_table)

        # Step 6: Full Contig Extraction
        print("[6/6] Extracting full original contigs for filtered eukaryotes...")
        extract_full_contigs(args.assembly, df, args.out_dir / "Contigs_with_18s.fasta")

        print(f"\nPipeline complete! All outputs saved to: {args.out_dir.resolve()}")

if __name__ == "__main__":
    main()
