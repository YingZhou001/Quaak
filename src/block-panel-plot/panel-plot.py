#!/usr/bin/env python3
"""
Plot contig mapping visualization from a whitespace-delimited .mat.block.csv file.

Input block format (no header), 10 columns:
  1  target_contig
  2  block_idx
  3  target_from
  4  target_to
  5  ref_contig
  6  ref_from
  7  ref_to
  8  pos_kmer
  9  neg_kmer
  10 boundary_kmers

Strand is derived from pos_kmer vs neg_kmer.

Usage:
  python panel-plot.py <input> [output.pdf] [--highlight FILE]
                       [--compress FLOAT] [--margin-frac FLOAT]

Coordinate scaling:
  Bar coordinates pass through a piecewise-linear scale map that keeps regions
  within `--margin-frac * span` of each block endpoint at 1:1 scale and
  compresses the rest by `--compress`. Use --compress 1 for plain linear
  layout (no compression); the default compresses empty gaps to make small
  blocks more visible.
"""

import sys
import re
import math
import os
import gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages


# ── colours ────────────────────────────────────────────────────────────────
FWD_COLOR  = '#3a7ebf'   # blue  – forward strand (pos_kmer > neg_kmer)
REV_COLOR  = '#bf3a3a'   # red   – reverse strand (neg_kmer > pos_kmer)
REF_BAR    = '#c8e6fa'   # light blue background bar (reference / target)
QRY_BAR    = '#faefc8'   # light yellow background bar (query)


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def parse_block_file(filename: str) -> pd.DataFrame:
    """Parse whitespace-separated .mat.block.csv (10 columns, no header).

    Columns: target_contig block_idx target_from target_to
             ref_contig ref_from ref_to pos_kmer neg_kmer boundary_kmers
    """
    sample_id = re.sub(r'\.mat\.block\.csv(\.gz)?$', '', os.path.basename(filename))
    rows = []
    with _open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) < 9:
                continue
            rows.append({
                'sample_id': sample_id,
                'query_ctg': parts[0],
                'idx':       int(parts[1]),
                'q_start':   int(parts[2]),
                'q_end':     int(parts[3]),
                't_ctg':     parts[4],
                't_start':   int(parts[5]),
                't_end':     int(parts[6]),
                'n_pos':     int(parts[7]),
                'n_neg':     int(parts[8]),
            })
    return pd.DataFrame(rows)


def choose_unit(values) -> tuple:
    """
    Pick a single unit (Mb, kb, or bp) for a collection of positions,
    based on the largest value, and return (divisor, suffix).
    """
    vmax = max(abs(v) for v in values if v is not None)
    if vmax >= 1_000_000:
        return 1_000_000, 'Mb'
    if vmax >= 1_000:
        return 1_000, 'kb'
    return 1, 'bp'


def fmt_with_unit(x: float, divisor: float, suffix: str) -> str:
    val = x / divisor
    if divisor == 1:
        return f'{int(val)} {suffix}'
    return f'{val:.2f} {suffix}'


def build_scale_map(intervals, total_min, total_max,
                    margin_frac=0.005, compress=0.001):
    """Piecewise-linear map from data coords to [0, 1] that compresses regions
    far from any block boundary. Regions within `margin_frac * span` of a
    block endpoint keep scale 1; other regions are scaled by `compress`.
    With compress=1 the mapping is exactly linear.
    """
    span   = max(total_max - total_min, 1)
    margin = margin_frac * span

    breakpoints = sorted(set(
        [total_min, total_max] +
        [x for lo, hi in intervals for x in (lo, hi)]
    ))

    zone_xs = {total_min, total_max}
    for bp in breakpoints:
        zone_xs.add(max(bp - margin, total_min))
        zone_xs.add(min(bp + margin, total_max))
    zone_xs = sorted(zone_xs)

    key_data = [total_min]
    key_cum  = [0.0]
    cumulative = 0.0
    for x0, x1 in zip(zone_xs, zone_xs[1:]):
        mid   = (x0 + x1) / 2
        near  = any(abs(mid - bp) <= margin for bp in breakpoints)
        scale = 1.0 if near else compress
        cumulative += (x1 - x0) * scale
        key_data.append(x1)
        key_cum.append(cumulative)

    total    = key_cum[-1] or 1.0
    key_norm = [v / total for v in key_cum]
    return list(zip(key_data, key_norm))


def apply_scale_map(x, scale_map):
    xs = [p[0] for p in scale_map]
    ys = [p[1] for p in scale_map]
    return float(np.interp(x, xs, ys))


def draw_subplot(ax, blocks: pd.DataFrame,
                 q_range: tuple, t_range: tuple, title: str,
                 compress: float = 0.001,
                 margin_frac: float = 0.005) -> None:
    """
    Draw one mapping panel on *ax*.

    Layout (normalised x ∈ [0,1], y ∈ [0,1]):
      • Reference bar  at top    (y ~ 0.75–0.83)
      • Query bar      at bottom (y ~ 0.15–0.23)
      • Trapezoid ribbons connect matching blocks between the two bars.

    Strand: n_pos vs n_neg (forward ribbons parallel, reverse ribbons cross).

    Coordinates pass through build_scale_map — gaps between blocks are
    compressed by `compress` (set compress=1 for plain linear scaling).
    """
    q_min, q_max = q_range
    t_min, t_max = t_range
    q_span = max(q_max - q_min, 1)
    t_span = max(t_max - t_min, 1)

    q_intervals = list(zip(blocks['q_start'], blocks['q_end']))
    t_intervals = list(zip(blocks['t_start'], blocks['t_end']))
    q_map = build_scale_map(q_intervals, q_min, q_max, margin_frac, compress)
    t_map = build_scale_map(t_intervals, t_min, t_max, margin_frac, compress)

    # Use the larger span to set relative bar widths so the visual length
    # difference between ref and query is preserved.
    total_span = max(q_span, t_span)
    t_width = t_span / total_span
    q_width = q_span / total_span

    def nq(x):
        return apply_scale_map(x, q_map) * q_width

    def nt(x):
        return apply_scale_map(x, t_map) * t_width

    REF_BOT = 0.75
    QRY_BOT = 0.15
    BAR_H   = 0.08
    REF_TOP = REF_BOT + BAR_H
    QRY_TOP = QRY_BOT + BAR_H

    # Background bars (each drawn only as wide as its own span)
    ax.add_patch(patches.Rectangle((0, REF_BOT), t_width, BAR_H,
                                   color=REF_BAR, zorder=2))
    ax.add_patch(patches.Rectangle((0, QRY_BOT), q_width, BAR_H,
                                   color=QRY_BAR, zorder=2))

    for _, blk in blocks.iterrows():
        n_pos = blk['n_pos']
        n_neg = blk['n_neg']

        if n_pos > n_neg:
            color, strand = FWD_COLOR, '+'
        elif n_neg > n_pos:
            color, strand = REV_COLOR, '-'
        else:
            continue   # ambiguous – skip
        ribbon_alpha, bar_alpha = 0.30, 0.85

        ts = nt(blk['t_start'])
        te = nt(blk['t_end'])
        qs = nq(blk['q_start'])
        qe = nq(blk['q_end'])

        # Ribbon polygon:
        #   forward:  ref_left→query_left , ref_right→query_right  (parallel)
        #   reverse:  ref_left→query_right, ref_right→query_left   (crossing)
        if strand == '+':
            poly_pts = [(ts, REF_BOT), (te, REF_BOT),
                        (qe, QRY_TOP), (qs, QRY_TOP)]
        else:
            poly_pts = [(ts, REF_BOT), (te, REF_BOT),
                        (qs, QRY_TOP), (qe, QRY_TOP)]

        ax.add_patch(plt.Polygon(poly_pts, closed=True,
                                 facecolor=color, edgecolor=color,
                                 alpha=ribbon_alpha, linewidth=0.4, zorder=1))

        # Highlight the block interval on each bar
        ax.add_patch(patches.Rectangle((ts, REF_BOT), te - ts, BAR_H,
                                       color=color, alpha=bar_alpha, zorder=3))
        ax.add_patch(patches.Rectangle((qs, QRY_BOT), qe - qs, BAR_H,
                                       color=color, alpha=bar_alpha, zorder=3))

    # Lane labels (left of bars)
    ax.text(-0.02, REF_BOT + BAR_H / 2, 'Ref',
            ha='right', va='center', fontsize=10, color='#333')
    ax.text(-0.02, QRY_BOT + BAR_H / 2, 'Query',
            ha='right', va='center', fontsize=10, color='#333')

    # Coordinate tick labels – single unit across both lanes
    div, suf = choose_unit([q_min, q_max, t_min, t_max])
    for x, val in [(0, t_min), (t_width, t_max)]:
        ax.text(x, REF_BOT - 0.03, fmt_with_unit(val, div, suf),
                ha='left' if x == 0 else 'right',
                va='top', fontsize=8, color='#555')
    for x, val in [(0, q_min), (q_width, q_max)]:
        ax.text(x, QRY_BOT - 0.03, fmt_with_unit(val, div, suf),
                ha='left' if x == 0 else 'right',
                va='top', fontsize=8, color='#555')

    ax.set_xlim(-0.14, 1.08)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title(title, fontsize=9, pad=2)


def main():
    import argparse
    ap = argparse.ArgumentParser(description='Plot contig mapping visualization.')
    ap.add_argument('input',  nargs='?', default='difficult_ctg.txt')
    ap.add_argument('output', nargs='?', default='mapping_plot.pdf')
    ap.add_argument('--highlight', default='manual-rescue.csv',
                    help='Optional CSV whose first two columns are query_ctg and t_ctg to highlight (default: manual-rescue.csv)')
    ap.add_argument('--compress', type=float, default=0.001,
                    help='Scale factor applied to regions between blocks; '
                         '1.0 disables compression for plain linear scaling (default: 0.001)')
    ap.add_argument('--margin-frac', type=float, default=0.005,
                    help='Fraction of total span kept at 1:1 scale around each block endpoint (default: 0.005)')
    args = ap.parse_args()
    infile  = args.input
    outfile = args.output

    rescue_pairs = {}  # (query_ctg, t_ctg) → box color
    if os.path.exists(args.highlight):
        with open(args.highlight) as fh:
            for line in fh:
                parts = line.strip().split()
                if len(parts) >= 3:
                    strand = 1 if float(parts[2]) > 0 else -1
                    color  = FWD_COLOR if strand == 1 else REV_COLOR
                    rescue_pairs[(parts[0], parts[1])] = color
        print(f"Loaded {len(rescue_pairs)} highlight pairs from {args.highlight}.")

    df = parse_block_file(infile)
    print(f"Parsed {len(df)} rows (block format).")
    blocks_df = df.copy()

    # One subplot per (query_ctg, t_ctg) pair.
    # Preserve input-file order: key on first-appearance index of each pair
    pair_order = (
        blocks_df[['query_ctg', 't_ctg']]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    groups = [
        ((row.query_ctg, row.t_ctg),
         blocks_df[(blocks_df['query_ctg'] == row.query_ctg) &
                   (blocks_df['t_ctg']     == row.t_ctg)])
        for _, row in pair_order.iterrows()
    ]

    def _chr_sort_key(item):
        """Natural sort on t_ctg (chr2 before chr10), then by t_start."""
        t_ctg = item[0][1]
        t_pos = item[1]['t_start'].min()
        natural = tuple(int(p) if p.isdigit() else p
                        for p in re.split(r'(\d+)', t_ctg))
        return (natural, t_pos)

    groups.sort(key=_chr_sort_key)
    n_groups = len(groups)
    print(f"Found {n_groups} unique (query_ctg, target_ctg) pairs after filtering.")

    COLS     = 4
    ROWS     = 3
    PER_PAGE = COLS * ROWS
    n_pages  = math.ceil(n_groups / PER_PAGE)

    # Legend patches
    fwd_patch = patches.Patch(color=FWD_COLOR, alpha=0.8, label='Forward')
    rev_patch = patches.Patch(color=REV_COLOR, alpha=0.8, label='Reverse')

    with PdfPages(outfile) as pdf:
        for page in range(n_pages):
            fig, axes = plt.subplots(ROWS, COLS, figsize=(18, 10))
            axes_flat = axes.flatten()

            for slot in range(PER_PAGE):
                gi  = page * PER_PAGE + slot
                ax  = axes_flat[slot]
                if gi >= n_groups:
                    ax.axis('off')
                    continue

                (query_ctg, t_ctg), grp = groups[gi]

                # Range from actual blocks (covers all plottable and uncertain blocks)
                q_min = grp['q_start'].min()
                q_max = grp['q_end'].max()
                t_min = grp['t_start'].min()
                t_max = grp['t_end'].max()

                title = f"{query_ctg} {t_ctg}"

                draw_subplot(ax, grp, (q_min, q_max), (t_min, t_max), title,
                             compress=args.compress,
                             margin_frac=args.margin_frac)

                if (query_ctg, t_ctg) in rescue_pairs:
                    box_color = rescue_pairs[(query_ctg, t_ctg)]
                    ax.add_patch(patches.FancyBboxPatch(
                        (0, 0), 1, 0.92,
                        transform=ax.transAxes,
                        boxstyle='square,pad=0.02',
                        fill=False, edgecolor=box_color,
                        linewidth=2.5, linestyle='--',
                        zorder=10, clip_on=False,
                    ))

            # Shared legend on each page
            fig.legend(handles=[fwd_patch, rev_patch],
                       loc='lower center', ncol=2, fontsize=11,
                       frameon=False, bbox_to_anchor=(0.5, 0.0))

            fig.suptitle(f'Contig Block Mapping  (page {page + 1}/{n_pages})',
                         fontsize=13, y=1.01)
            plt.tight_layout(rect=[0, 0.04, 1, 1])
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            print(f"  Page {page + 1}/{n_pages} written.")

    print(f"\nDone. Output saved to: {outfile}")


if __name__ == '__main__':
    main()
