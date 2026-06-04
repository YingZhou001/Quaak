#!/usr/bin/env python3
"""
Plot anchor-flanked region visualization from sel.annot.tsv.

Input format (tab-separated, no header), 10 columns:
  1  query_ctg
  2  q_start            (NA for DEL rows)
  3  q_end              (NA for DEL rows)
  4  t_ctg
  5  t_start
  6  t_end
  7  t_outer_start
  8  t_outer_end
  9  count
  10 type               ANC | INV | DUP | DEL | other

Each subplot shows one "anchor + (any) + anchor" region:
two consecutive ANC rows (within the same query/target contig pair)
together with all non-ANC rows between them.

Usage:
  python anchor-plot.py <input> [output.pdf]
                        [--compress FLOAT] [--margin-frac FLOAT]
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
FWD_COLOR  = '#3a7ebf'   # blue   – ANC (forward anchor)
REV_COLOR  = '#bf3a3a'   # red    – INV (inversion)
DUP_COLOR  = '#7a3abf'   # purple – DUP
DEL_COLOR  = '#555555'   # grey   – DEL (ref-only)
OTHER_COLOR = '#bf8a3a'  # orange – other
REF_BAR    = '#c8e6fa'   # light blue background bar (reference)
QRY_BAR    = '#faefc8'   # light yellow background bar (query)

TYPE_COLOR = {
    'ANC':   FWD_COLOR,
    'INV':   REV_COLOR,
    'DUP':   DUP_COLOR,
    'DEL':   DEL_COLOR,
    'other': OTHER_COLOR,
}


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def _maybe_int(x):
    if x is None or x == '' or x.upper() == 'NA':
        return None
    return int(x)


def parse_annot_file(filename: str) -> pd.DataFrame:
    """Parse the 10-column sel.annot.tsv file (no header)."""
    rows = []
    with _open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 10:
                parts = line.strip().split()
            if len(parts) < 10:
                continue
            rows.append({
                'query_ctg':     parts[0],
                'q_start':       _maybe_int(parts[1]),
                'q_end':         _maybe_int(parts[2]),
                't_ctg':         parts[3],
                't_start':       _maybe_int(parts[4]),
                't_end':         _maybe_int(parts[5]),
                't_outer_start': _maybe_int(parts[6]),
                't_outer_end':   _maybe_int(parts[7]),
                'count':         _maybe_int(parts[8]) or 0,
                'type':          parts[9],
            })
    return pd.DataFrame(rows)


def choose_unit(values) -> tuple:
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
    """Piecewise-linear map from data coords to [0, 1] that compresses
    regions far from any block boundary. Regions within
    `margin_frac * span` of an endpoint keep scale 1; other regions are
    scaled by `compress`. With compress=1 the mapping is exactly linear.
    """
    span   = max(total_max - total_min, 1)
    margin = margin_frac * span

    breakpoints = sorted(set(
        [total_min, total_max] +
        [x for lo, hi in intervals for x in (lo, hi) if x is not None]
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
                 q_range: tuple, t_range: tuple,
                 compress: float = 0.001,
                 margin_frac: float = 0.005) -> None:
    """Draw one anchor-flanked region on *ax*."""
    q_min, q_max = q_range
    t_min, t_max = t_range
    q_span = max(q_max - q_min, 1)
    t_span = max(t_max - t_min, 1)

    q_intervals = [(s, e) for s, e in zip(blocks['q_start'], blocks['q_end'])
                   if s is not None and e is not None]
    t_intervals = [(s, e) for s, e in zip(blocks['t_start'], blocks['t_end'])
                   if s is not None and e is not None]
    q_map = build_scale_map(q_intervals, q_min, q_max, margin_frac, compress)
    t_map = build_scale_map(t_intervals, t_min, t_max, margin_frac, compress)

    total_span = max(q_span, t_span)
    t_width = t_span / total_span
    q_width = q_span / total_span

    # Detect crossed flanking anchors: if the first and last ANC rows run in
    # opposite order on the query vs. the ref (e.g. q[1,5]->t[21,25] together
    # with q[20,24]->t[1,5]), the query segment is reverse-oriented. Mirror the
    # query axis so the anchors read parallel and true inversions stay crossed.
    def _is_anc(t):
        return 'ANC' in [s.strip() for s in str(t).split(',')]

    q_reversed = False
    anc = blocks[blocks['type'].apply(_is_anc)]
    if len(anc) >= 2:
        a0, a1 = anc.iloc[0], anc.iloc[-1]
        coords = [a0['q_start'], a0['q_end'], a1['q_start'], a1['q_end'],
                  a0['t_start'], a0['t_end'], a1['t_start'], a1['t_end']]
        if all(c is not None for c in coords):
            dq = (a1['q_start'] + a1['q_end']) - (a0['q_start'] + a0['q_end'])
            dt = (a1['t_start'] + a1['t_end']) - (a0['t_start'] + a0['t_end'])
            q_reversed = dq * dt < 0

    def nq(x):
        v = apply_scale_map(x, q_map)
        if q_reversed:
            v = 1.0 - v
        return v * q_width

    def nt(x):
        return apply_scale_map(x, t_map) * t_width

    REF_BOT = 0.75
    QRY_BOT = 0.15
    BAR_H   = 0.08
    REF_TOP = REF_BOT + BAR_H
    QRY_TOP = QRY_BOT + BAR_H

    ax.add_patch(patches.Rectangle((0, REF_BOT), t_width, BAR_H,
                                   color=REF_BAR, zorder=2))
    ax.add_patch(patches.Rectangle((0, QRY_BOT), q_width, BAR_H,
                                   color=QRY_BAR, zorder=2))

    for _, blk in blocks.iterrows():
        types = [t.strip() for t in str(blk['type']).split(',') if t.strip()]
        is_inv = 'INV' in types
        color = REV_COLOR if is_inv else FWD_COLOR
        ribbon_alpha, bar_alpha = 0.30, 0.85

        ts = nt(blk['t_start']) if blk['t_start'] is not None else None
        te = nt(blk['t_end'])   if blk['t_end']   is not None else None
        qs = nq(blk['q_start']) if blk['q_start'] is not None else None
        qe = nq(blk['q_end'])   if blk['q_end']   is not None else None

        # DEL: only ref-side coordinates available, no ribbon
        if qs is None or qe is None:
            if ts is not None and te is not None:
                ax.add_patch(patches.Rectangle((ts, REF_BOT), te - ts, BAR_H,
                                               color=color, alpha=bar_alpha,
                                               zorder=3))
            continue

        # Ribbon orientation: an INV block crosses; others are parallel. When
        # the query axis is mirrored the vertex handedness flips, so swap the
        # shape (is_inv XOR q_reversed) to keep the visible crossing == is_inv.
        if is_inv ^ q_reversed:
            poly_pts = [(ts, REF_BOT), (te, REF_BOT),
                        (qs, QRY_TOP), (qe, QRY_TOP)]
        else:
            poly_pts = [(ts, REF_BOT), (te, REF_BOT),
                        (qe, QRY_TOP), (qs, QRY_TOP)]

        ax.add_patch(plt.Polygon(poly_pts, closed=True,
                                 facecolor=color, edgecolor=color,
                                 alpha=ribbon_alpha, linewidth=0.4, zorder=1))

        ax.add_patch(patches.Rectangle((ts, REF_BOT), te - ts, BAR_H,
                                       color=color, alpha=bar_alpha, zorder=3))
        ax.add_patch(patches.Rectangle((min(qs, qe), QRY_BOT), abs(qe - qs),
                                       BAR_H, color=color, alpha=bar_alpha,
                                       zorder=3))

    ax.text(-0.02, REF_BOT + BAR_H / 2, 'Ref',
            ha='right', va='center', fontsize=16, color='#333')
    ax.text(-0.02, QRY_BOT + BAR_H / 2, 'Query',
            ha='right', va='center', fontsize=16, color='#333')

    div, suf = choose_unit([q_min, q_max, t_min, t_max])
    for x, val in [(0, t_min), (t_width, t_max)]:
        ax.text(x, REF_BOT - 0.03, fmt_with_unit(val, div, suf),
                ha='left' if x == 0 else 'right',
                va='top', fontsize=14, color='#555')
    q_ends = [(0, q_max), (q_width, q_min)] if q_reversed \
        else [(0, q_min), (q_width, q_max)]
    for x, val in q_ends:
        ax.text(x, QRY_BOT - 0.03, fmt_with_unit(val, div, suf),
                ha='left' if x == 0 else 'right',
                va='top', fontsize=14, color='#555')

    ax.set_xlim(-0.14, 1.08)
    ax.set_ylim(0, 1)
    ax.axis('off')


def build_anchor_groups(df: pd.DataFrame):
    """Return list of DataFrames, each a region of interest.

    Anchors are rows with type == 'ANC'. Within each (query_ctg, t_ctg):
      * every consecutive pair of anchors yields anchor + (any) + anchor;
      * non-ANC blocks before the first anchor or after the last anchor
        (or around a lone anchor) yield a one-sided region whose open
        boundary is taken from the block span itself;
      * a contig with blocks but no anchor yields that block span.
    """
    groups = []
    for (qctg, tctg), sub in df.groupby(['query_ctg', 't_ctg'], sort=False):
        sub = sub.reset_index(drop=True)
        anchor_idx = sub.index[sub['type'] == 'ANC'].tolist()

        if not anchor_idx:
            if len(sub):
                groups.append(sub.copy())
            continue

        # blocks before the first anchor (missing left anchor)
        first = anchor_idx[0]
        if first > 0:
            groups.append(sub.iloc[:first + 1].copy())

        # anchor + (any) + anchor for each consecutive pair
        for a, b in zip(anchor_idx, anchor_idx[1:]):
            groups.append(sub.iloc[a:b + 1].copy())

        # blocks after the last anchor (missing right anchor)
        last = anchor_idx[-1]
        if last < len(sub) - 1:
            groups.append(sub.iloc[last:].copy())

    return groups


def anchor_gap(group: pd.DataFrame):
    """Return (start, end) of the ref gap *between* the two flanking anchors.

    The gap is the space between the inner edges of the left and right ANC
    blocks: left anchor's far edge -> right anchor's near edge. Returns
    (None, None) when the group has fewer than two anchors with ref coords.
    """
    anc = group[group['type'].str.contains('ANC', na=False)]
    ivals = [(min(s, e), max(s, e))
             for s, e in zip(anc['t_start'], anc['t_end'])
             if s is not None and e is not None]
    if len(ivals) < 2:
        return None, None
    left, right = min(ivals), max(ivals)
    return left[1], right[0]


def overlay_dels(group: pd.DataFrame, dels: pd.DataFrame) -> pd.DataFrame:
    """Add any DEL rows whose ref interval overlaps this group's ref range."""
    if dels.empty:
        return group
    qctg = group['query_ctg'].iloc[0]
    tctg = group['t_ctg'].iloc[0]
    t_lo = group['t_start'].min()
    t_hi = group['t_end'].max()
    cand = dels[(dels['query_ctg'] == qctg) & (dels['t_ctg'] == tctg)]
    overlap = cand[(cand['t_end']   >= t_lo) &
                   (cand['t_start'] <= t_hi)]
    if overlap.empty:
        return group
    return pd.concat([group, overlap], ignore_index=True)


def choose_grid(n, rows=None, cols=None, max_cols=5):
    """Pick a panel grid for *n* regions.

    Honours whichever of *rows*/*cols* the user fixed and derives the rest.
    With neither set, returns a near-square grid (cols capped at *max_cols*),
    so a single page is used when the regions fit and paging kicks in beyond.
    """
    n = max(n, 1)
    if rows and cols:
        return rows, cols
    if cols:
        return math.ceil(n / cols), cols
    if rows:
        return rows, math.ceil(n / rows)
    cols = min(max_cols, max(1, math.ceil(math.sqrt(n))))
    rows = math.ceil(n / cols)
    return rows, cols


def main():
    import argparse
    ap = argparse.ArgumentParser(
        description='Plot anchor-flanked regions from sel.annot.tsv.')
    ap.add_argument('input',  nargs='?', default='sel.annot.tsv')
    ap.add_argument('output', nargs='?', default='anchor_plot.pdf')
    ap.add_argument('--compress', type=float, default=0.001,
                    help='Scale factor applied to regions between blocks; '
                         '1.0 disables compression (default: 0.001)')
    ap.add_argument('--margin-frac', type=float, default=0.005,
                    help='Fraction of total span kept at 1:1 scale around '
                         'each block endpoint (default: 0.005)')
    ap.add_argument('--rows', type=int, default=None,
                    help='Panels per page, vertically '
                         '(default: auto from region count)')
    ap.add_argument('--cols', type=int, default=None,
                    help='Panels per page, horizontally '
                         '(default: auto from region count)')
    ap.add_argument('--max-cols', type=int, default=5,
                    help='Upper bound on columns when auto-selecting the grid '
                         '(default: 5)')
    ap.add_argument('--panel-w', type=float, default=6.0,
                    help='Width of each panel in inches (default: 6.0)')
    ap.add_argument('--panel-h', type=float, default=3.0,
                    help='Height of each panel in inches (default: 3.0)')
    args = ap.parse_args()

    df = parse_annot_file(args.input)
    print(f"Parsed {len(df)} rows from {args.input}.")

    df = df[df['type'] != 'DEL'].copy()
    groups = build_anchor_groups(df)
    n_groups = len(groups)
    print(f"Built {n_groups} anchor+...+anchor regions.")
    if n_groups == 0:
        print("Nothing to plot.")
        return

    nrows, ncols = choose_grid(n_groups, args.rows, args.cols, args.max_cols)
    panel_w, panel_h = args.panel_w, args.panel_h
    per_page = nrows * ncols
    n_pages = math.ceil(n_groups / per_page)

    with PdfPages(args.output) as pdf:
        for page in range(n_pages):
            page_groups = groups[page * per_page:(page + 1) * per_page]
            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(panel_w * ncols, panel_h * nrows))
            axes_flat = np.array(axes).flatten()
            for ax in axes_flat[len(page_groups):]:
                ax.axis('off')

            for grp, ax in zip(page_groups, axes_flat):
                q_vals = [v for v in list(grp['q_start']) + list(grp['q_end'])
                          if v is not None]
                t_vals = [v for v in list(grp['t_start']) + list(grp['t_end'])
                          if v is not None]
                q_min, q_max = (min(q_vals), max(q_vals)) if q_vals else (0, 1)
                t_min, t_max = (min(t_vals), max(t_vals)) if t_vals else (0, 1)

                draw_subplot(ax, grp, (q_min, q_max), (t_min, t_max),
                             compress=args.compress,
                             margin_frac=args.margin_frac)

                qctg = grp['query_ctg'].iloc[0]
                tctg = grp['t_ctg'].iloc[0]
                gap_start, gap_end = anchor_gap(grp)
                if gap_start is None:
                    span = f'{t_min} - {t_max}'
                else:
                    span = f'{gap_start} - {gap_end}'
                ax.set_title(f'{qctg} {tctg}:{span}',
                             fontsize=13, color='#222')

            plt.tight_layout(rect=[0, 0.02, 1, 1])
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
        print(f"  {n_pages} page(s) with {n_groups} panels written "
              f"({nrows}x{ncols} per page).")

    print(f"\nDone. Output saved to: {args.output}")


if __name__ == '__main__':
    main()
