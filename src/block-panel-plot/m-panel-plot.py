#!/usr/bin/env python3
"""
Multi-mode block visualization.

Input: same 10-column whitespace-delimited block file as panel-plot.py
       (target_contig block_idx target_from target_to ref_contig
        ref_from ref_to pos_kmer neg_kmer boundary_kmers).
       Column 1 is treated as the contig (query); column 5 as the reference.

Modes (mutually exclusive, exactly one required):
  --ref-centric   one panel per reference; each panel shows
                    top    : the reference bar
                    middle : every contig that has at least one block on this
                             reference (so a contig mapping to N refs appears
                             on N panels), laid out at this-ref coords with
                             sweep-line lane packing for overlapping contigs.
                             Contigs whose MAJOR ref is this panel ("primary")
                             get a full-alpha backbone; contigs that also map
                             elsewhere ("shared") get a dimmed backbone.
                             Off-panel mappings extend the backbone length so
                             the contig is drawn at its true length, but no
                             secondary blocks/ribbons are drawn on this panel.

  --ctg-centric   one panel per contig; each panel shows
                    bottom : the contig bar (q-coords)
                    middle : the contig's major reference, blocks at q-coords
                    top    : blocks mapping to any secondary reference,
                             lane-packed by overlap on q-coords (omitted
                             if no secondary mappings exist)

Major reference per contig = the t_ctg with the largest total (n_pos + n_neg)
across all blocks of that contig; ties broken by first appearance order.
Use --min-kmer to drop blocks whose individual kmer count is below a
threshold (default 0, i.e. keep all).
"""

import argparse
import gzip
import math
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages


FWD_COLOR = '#3a7ebf'
REV_COLOR = '#e07b00'
REF_BAR   = '#c8e6fa'
CTG_BAR   = '#faefc8'
SEC_BAR   = '#e6d8f5'

# UCSC Giemsa stain colors for cytobands
CYTOBAND_COLORS = {
    'gneg':    '#ffffff',
    'gpos25':  '#cccccc',
    'gpos50':  '#999999',
    'gpos75':  '#666666',
    'gpos100': '#000000',
    'acen':    '#b03a3a',
    'gvar':    '#dcdcec',
    'stalk':   '#c8b6e0',
}

# Vertical fraction of a band/lane actually occupied by the bar
# (the rest is whitespace above/below). Reduce to make bars thinner.
BAR_FILL  = 0.35



def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def parse_strand_tsv(filename: str) -> dict:
    """Parse a strand file with columns: col1=ctg, col2=ref, col3=tau.
       Returns {(ctg, ref): tau}. Pairs with tau >= 0 stay forward; tau < 0
       means the contig should be reverse-complemented for that ref."""
    out = {}
    with _open(filename) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            p = line.split()
            if len(p) < 3:
                continue
            try:
                tau = float(p[2])
            except ValueError:
                continue
            out[(p[0], p[1])] = tau
    return out


def apply_strand_flips(df: pd.DataFrame, major: dict, strand: dict):
    """For each contig whose (ctg, major_ref) tau is < 0, reverse the contig:
    q_new = L - q_old (L = max q_end over the contig's blocks) and swap
    n_pos ↔ n_neg. All blocks of the contig flip together so the contig's
    q-frame stays consistent across major and secondary refs. Pairs missing
    from the strand file are left alone."""
    if not strand:
        return df, set()
    flipped = set()
    df = df.copy()
    for q, ref in major.items():
        tau = strand.get((q, ref))
        if tau is None or tau >= 0:
            continue
        flipped.add(q)
    if not flipped:
        return df, flipped
    for q in flipped:
        mask = df['query_ctg'] == q
        L = int(df.loc[mask, 'q_end'].max())
        old_qs = df.loc[mask, 'q_start'].copy()
        old_qe = df.loc[mask, 'q_end'].copy()
        df.loc[mask, 'q_start'] = L - old_qe
        df.loc[mask, 'q_end']   = L - old_qs
        old_pos = df.loc[mask, 'n_pos'].copy()
        df.loc[mask, 'n_pos'] = df.loc[mask, 'n_neg']
        df.loc[mask, 'n_neg'] = old_pos
    return df, flipped


def parse_cytobands(filename: str) -> dict:
    """Parse a UCSC cytoband BED → {chrom: [(start, end, name, stain), ...]}.
    Each chromosome is also indexed under the trailing 'chrN' token
    (e.g. 'CHM13#0#chr1' → 'chr1') so block files using bare chrom names
    still match."""
    bands = {}
    with _open(filename) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            p = line.rstrip().split('\t')
            if len(p) < 5:
                p = line.split()
            if len(p) < 5:
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                continue
            bands.setdefault(p[0], []).append((s, e, p[3], p[4]))
    out = {c: sorted(b) for c, b in bands.items()}
    for chrom in list(out.keys()):
        m = re.search(r'(chr[\dXYMxym]+)$', chrom)
        if m and m.group(1) not in out:
            out[m.group(1)] = out[chrom]
    return out


def lookup_cytobands(cytobands: dict, ref_name: str):
    if not cytobands:
        return None
    if ref_name in cytobands:
        return cytobands[ref_name]
    m = re.search(r'(chr[\dXYMxym]+)$', ref_name)
    if m and m.group(1) in cytobands:
        return cytobands[m.group(1)]
    return None


def centromere_span(bands):
    """(start, end) covering all 'acen' bands in `bands`, or None."""
    if not bands:
        return None
    s = e = None
    for bs, be, _, stain in bands:
        if stain != 'acen':
            continue
        s = bs if s is None else min(s, bs)
        e = be if e is None else max(e, be)
    return (s, e) if s is not None else None


def bands_span(bands):
    """(start, end) covering all bands in `bands`, or None."""
    if not bands:
        return None
    return min(b[0] for b in bands), max(b[1] for b in bands)


def parse_blocks(filename: str) -> pd.DataFrame:
    rows = []
    with _open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            p = line.strip().split()
            if len(p) < 9:
                continue
            rows.append({
                'query_ctg': p[0],
                'idx':       int(p[1]),
                'q_start':   int(p[2]),
                'q_end':     int(p[3]),
                't_ctg':     p[4],
                't_start':   int(p[5]),
                't_end':     int(p[6]),
                'n_pos':     int(p[7]),
                'n_neg':     int(p[8]),
            })
    return pd.DataFrame(rows)


def assign_major_ref(df: pd.DataFrame, min_kmer: int = 0,
                     min_pair_kmer: int = 0):
    """Return (df_filtered, {query_ctg: major_t_ctg}).

    `min_kmer`      drops blocks whose own (n_pos + n_neg) is below threshold.
    `min_pair_kmer` drops every block of a (query_ctg, t_ctg) pair whose
                    SUMMED kmer count across the pair is below threshold.
                    A contig whose major pair is dropped disappears entirely.
    """
    df = df.copy()
    df['kmer'] = df['n_pos'] + df['n_neg']
    if min_kmer > 0:
        df = df[df['kmer'] >= min_kmer].copy()

    first_seen = {}
    for i, (q, t) in enumerate(zip(df['query_ctg'], df['t_ctg'])):
        first_seen.setdefault((q, t), i)

    pair = (df.groupby(['query_ctg', 't_ctg'], as_index=False)['kmer']
              .sum())

    if min_pair_kmer > 0:
        keep = set((q, t) for q, t, k in
                   zip(pair['query_ctg'], pair['t_ctg'], pair['kmer'])
                   if k >= min_pair_kmer)
        df = df[df.apply(lambda r: (r['query_ctg'], r['t_ctg']) in keep,
                         axis=1)].copy()
        pair = pair[pair.apply(
            lambda r: (r['query_ctg'], r['t_ctg']) in keep, axis=1)]

    pair['order'] = [first_seen[(q, t)]
                     for q, t in zip(pair['query_ctg'], pair['t_ctg'])]
    pair = pair.sort_values(['query_ctg', 'kmer', 'order'],
                            ascending=[True, False, True])
    major = (pair.drop_duplicates('query_ctg', keep='first')
                 .set_index('query_ctg')['t_ctg']
                 .to_dict())
    return df, major


def strand_color(n_pos: int, n_neg: int):
    if n_pos > n_neg:
        return FWD_COLOR, '+'
    if n_neg > n_pos:
        return REV_COLOR, '-'
    return None, None


def lane_pack(intervals):
    """Sweep-line lane assignment, stable on input order for ties.

    intervals : list of (start, end, payload)
    returns   : (list of (start, end, payload, lane), n_lanes)

    An item shares a lane only if its start is strictly past the lane's
    current end (so [1,15] and [10,20] go to different lanes).
    """
    indexed = list(enumerate(intervals))
    indexed.sort(key=lambda x: (x[1][0], x[0]))

    lane_end = []
    out = [None] * len(intervals)
    for orig_i, (s, e, payload) in indexed:
        placed = False
        for li, end in enumerate(lane_end):
            if s >= end:
                lane_end[li] = e
                out[orig_i] = (s, e, payload, li)
                placed = True
                break
        if not placed:
            lane_end.append(e)
            out[orig_i] = (s, e, payload, len(lane_end) - 1)
    return out, max(len(lane_end), 1)


def choose_unit(values):
    vmax = max(abs(v) for v in values if v is not None)
    if vmax >= 1_000_000:
        return 1_000_000, 'Mb'
    if vmax >= 1_000:
        return 1_000, 'kb'
    return 1, 'bp'


def fmt_unit(x, divisor, suffix):
    val = x / divisor
    return f'{int(val)} {suffix}' if divisor == 1 else f'{val:.2f} {suffix}'


def natural_key(s):
    return tuple(int(p) if p.isdigit() else p
                 for p in re.split(r'(\d+)', str(s)))


# ── geometry helpers ──────────────────────────────────────────────────────

def _norm(x, lo, span):
    return (x - lo) / span if span > 0 else 0.0


def _ribbon(ax, x1a, x1b, y1, x2a, x2b, y2, color, strand):
    """Trapezoid: (x1a,y1)-(x1b,y1)-(x2b,y2)-(x2a,y2). Reverse strand swaps
    the bottom edge endpoints so the polygon crosses."""
    if strand == '-':
        pts = [(x1a, y1), (x1b, y1), (x2a, y2), (x2b, y2)]
    else:
        pts = [(x1a, y1), (x1b, y1), (x2b, y2), (x2a, y2)]
    ax.add_patch(plt.Polygon(pts, closed=True,
                             facecolor=color, edgecolor=color,
                             alpha=0.30, linewidth=0.4, zorder=1))


def _bar(ax, x, y, w, h, color, alpha=1.0, z=2):
    ax.add_patch(patches.Rectangle((x, y), w, h,
                                   color=color, alpha=alpha, zorder=z))


def _shrink(bot, top, fill=BAR_FILL):
    """Centered vertical sub-range covering `fill` fraction of [bot, top]."""
    pad = (top - bot) * (1 - fill) / 2
    return bot + pad, top - pad


def _lane_y(band_top, band_bot, n_lanes, li, fill=0.8):
    """Lane spans the full band height divided by n_lanes (top-down).
    `fill` is the bar's fraction of that lane (rest is padding)."""
    lane_h = (band_top - band_bot) / max(n_lanes, 1)
    top = band_top - li * lane_h
    bot = top - lane_h
    pad = lane_h * (1 - fill) / 2
    return bot + pad, top - pad


def _lane_y_up(band_bot, band_top, n_lanes, li, fill=0.8):
    """Same as _lane_y but lane 0 sits at the bottom (bottom-up stacking)."""
    lane_h = (band_top - band_bot) / max(n_lanes, 1)
    bot = band_bot + li * lane_h
    top = bot + lane_h
    pad = lane_h * (1 - fill) / 2
    return bot + pad, top - pad


# ── ref-centric panel ─────────────────────────────────────────────────────

def _fit_name(name, max_len=48):
    """Middle-truncate a name with an ellipsis if it's longer than max_len."""
    if len(name) <= max_len:
        return name
    keep = max_len - 1
    head = keep // 2
    tail = keep - head
    return name[:head] + '…' + name[-tail:]


def _name_fontsize(name, base):
    """Shrink the font when the name gets long; keep base for short names."""
    n = len(name)
    if n > 36:
        return max(base - 4, 6)
    if n > 24:
        return max(base - 2, 7)
    return base


def draw_ref_panel(ax, ref_ctg, ref_blocks_by_ctg, sec_blocks_by_ctg,
                   x_lo, x_hi, global_span=None, labeled_contigs=None,
                   font_scale=1.0, one_plot=False, centromere=None,
                   primary_of=None, cytobands=None):
    """
    ref_blocks_by_ctg : {query_ctg: DataFrame of major-ref blocks}.
    sec_blocks_by_ctg : {query_ctg: DataFrame of secondary-ref blocks};
                        used only to extend the contig backbone to its full
                        q-range (no separate secondary band is drawn).
    x_lo, x_hi        : major-reference coord range
    global_span       : if given, normalize coords against this span instead
                        of (x_hi - x_lo) so that row widths reflect contig
                        length relative to the longest panel.
    """
    panel_span = max(x_hi - x_lo, 1)
    span = global_span if global_span else panel_span
    width_frac = panel_span / span

    # ── prepare lane assignments (needed before we can size bands) ──
    # Each contig is anchored on the ref axis at offset = min(t_start) of its
    # major-ref blocks. Within the contig lane, we plot in q-coords shifted
    # by that offset, so the backbone reflects actual contig length
    # (q_max - q_min) rather than the ref footprint.
    contigs = list(ref_blocks_by_ctg.keys())
    contig_offset = {}     # q-origin position on ref axis
    contig_qmin = {}
    contig_qmax = {}
    contig_intervals = []
    for q in contigs:
        mb = ref_blocks_by_ctg[q]
        sb = sec_blocks_by_ctg.get(q)
        anchor = mb.loc[mb['t_start'].idxmin()]
        t_anchor = int(anchor['t_start'])
        q_anchor = int(anchor['q_start'])
        q_lo = int(mb['q_start'].min())
        q_hi = int(mb['q_end'].max())
        if sb is not None and len(sb):
            q_lo = min(q_lo, int(sb['q_start'].min()))
            q_hi = max(q_hi, int(sb['q_end'].max()))
        offset = max(0, t_anchor - (q_anchor - q_lo))
        contig_offset[q] = offset
        contig_qmin[q] = q_lo
        contig_qmax[q] = q_hi
        contig_intervals.append((offset, offset + (q_hi - q_lo), q))
    packed, n_mid_lanes = lane_pack(contig_intervals)
    lane_of = {p[2]: p[3] for p in packed}

    # ── band layout (units: PANEL_UNIT per lane, GAP_UNIT between bands) ──
    u, g = float(PANEL_UNIT), float(GAP_UNIT)
    # When the caller restricts labels (e.g. --one-page picks just the top
    # contig), draw them as left-side row labels next to 'Ref' and skip the
    # bottom label strip entirely.
    use_left_labels = labeled_contigs is not None
    lbl_units = 0 if use_left_labels else LABEL_BAND_UNITS
    # cytoband strip sits ABOVE the ref bar (same x-extent / width)
    has_cyt = bool(cytobands)
    cyt_units = u * 0.3 if has_cyt else 0.0
    cyt_gap = g * 0.0 if has_cyt else 0.0
    total = (cyt_units + cyt_gap + u + g + n_mid_lanes * u
             + (0 if use_left_labels else g + lbl_units))

    # top-down stacking inside ax ylim [0, 1]
    CYT_TOP = 1.0
    CYT_BOT = CYT_TOP - cyt_units / total
    REF_TOP = CYT_BOT - cyt_gap / total
    REF_BOT = REF_TOP - u / total
    MID_TOP = REF_BOT - g / total
    MID_BOT = MID_TOP - n_mid_lanes * u / total
    LBL_TOP = MID_BOT - g / total
    LBL_BOT = LBL_TOP - lbl_units / total

    # Reference bar (thinner than the band it sits in). If cytobands are
    # given, the ref bar spans only the cytoband range; otherwise it spans
    # the full panel.
    rb_bot, rb_top = _shrink(REF_BOT, REF_TOP)
    if cytobands:
        bsp = bands_span(cytobands)
        rb_x0 = _norm(bsp[0], x_lo, span)
        rb_x1 = _norm(bsp[1], x_lo, span)
    else:
        rb_x0, rb_x1 = 0.0, width_frac
    _bar(ax, rb_x0, rb_bot, rb_x1 - rb_x0, rb_top - rb_bot, REF_BAR, z=2)

    # Cytoband strip above the ref bar — full ideogram clipped to panel range
    if has_cyt:
        for s, e, _name, stain in cytobands:
            ss = max(s, x_lo)
            ee = min(e, x_hi)
            if ee <= ss:
                continue
            color = CYTOBAND_COLORS.get(stain, '#dddddd')
            x0 = _norm(ss, x_lo, span)
            x1 = _norm(ee, x_lo, span)
            _bar(ax, x0, CYT_BOT, x1 - x0, CYT_TOP - CYT_BOT,
                 color, alpha=1.0, z=2)
        # outline around the strip so empty/light bands stay visible —
        # spans only the cytoband range, matching the ref bar.
        ax.add_patch(patches.Rectangle(
            (rb_x0, CYT_BOT), rb_x1 - rb_x0, CYT_TOP - CYT_BOT,
            fill=False, edgecolor='#444', linewidth=0.4, zorder=3))

    def lane_y(li):
        return _lane_y(MID_TOP, MID_BOT, n_mid_lanes, li, fill=BAR_FILL)

    # ── contig lanes + ribbons to ref bar ──
    for q in contigs:
        b = ref_blocks_by_ctg[q]
        li = lane_of[q]
        y_bot, y_top = lane_y(li)
        offset = contig_offset[q]
        q_min = contig_qmin[q]
        backbone_lo = offset
        backbone_hi = offset + (contig_qmax[q] - q_min)
        # contig backbone bar reflects actual contig length (q_max - q_min)
        bb_x0 = _norm(backbone_lo, x_lo, span)
        bb_x1 = _norm(backbone_hi, x_lo, span)
        is_primary = True if primary_of is None else primary_of.get(q, True)
        _bar(ax, bb_x0, y_bot, bb_x1 - bb_x0,
             y_top - y_bot, CTG_BAR,
             alpha=0.9 if is_primary else 0.35, z=2)
        # strike-through line grouping all blocks of this contig.
        # Deterministic per-contig jitter within the lane so adjacent lanes
        # don't visually merge.
        rng = np.random.default_rng(abs(hash(q)) % (2**32))
        lane_h = y_top - y_bot
        jitter = (rng.random() - 0.5) * lane_h * 0.5
        y_strike = (y_bot + y_top) / 2 + jitter
        ax.plot([bb_x0, bb_x1], [y_strike, y_strike],
                color='black', linewidth=0.6, zorder=10)
        # block-by-block: highlight on contig lane is in q-coords (shifted by
        # offset); on the ref bar it is in t-coords. Ribbon between them
        # encodes any length difference (indels) and strand (cross = inv).
        for _, blk in b.iterrows():
            color, strand = strand_color(blk['n_pos'], blk['n_neg'])
            if color is None:
                continue
            ts = _norm(blk['t_start'],            x_lo, span)
            te = _norm(blk['t_end'],              x_lo, span)
            qs = _norm(blk['q_start'] - q_min + offset, x_lo, span)
            qe = _norm(blk['q_end']   - q_min + offset, x_lo, span)
            _bar(ax, qs, y_bot, qe - qs, y_top - y_bot,
                 color, alpha=0.85, z=3)
            _bar(ax, ts, rb_bot, te - ts, rb_top - rb_bot,
                 color, alpha=0.85, z=3)
            _ribbon(ax,
                    ts, te, rb_bot,
                    qs, qe, y_top,
                    color, strand)

    label_set = (set(contigs) if labeled_contigs is None
                 else set(labeled_contigs))
    if use_left_labels:
        # Left-side row labels (same style as 'Ref'), one per LANE at the
        # lane's y-center. A lane may pack several non-overlapping contigs;
        # the highest-kmer contig in the lane represents it. Long names
        # shrink and then truncate.
        lane_rep = {}   # lane index -> (kmer_sum, contig)
        for q in contigs:
            li = lane_of[q]
            k = int(ref_blocks_by_ctg[q]['kmer'].sum())
            if li not in lane_rep or k > lane_rep[li][0]:
                lane_rep[li] = (k, q)
        for li, (_k, q) in lane_rep.items():
            y_bot, y_top = lane_y(li)
            shown = _fit_name(q, max_len=48)
            ax.text(-0.02, (y_bot + y_top) / 2, shown,
                    ha='right', va='center',
                    fontsize=_name_fontsize(q, 11), color='#333')
    else:
        # Bottom strip — names spaced evenly across the width, each connected
        # by a leader line up to its start position on the backbone.
        qs_sorted = [q for q in sorted(contigs, key=lambda q: contig_offset[q])
                     if q in label_set]
        n = len(qs_sorted)
        for i, q in enumerate(qs_sorted):
            slot_x = (i + 0.5) / n * width_frac if n else 0.0
            anchor_x = _norm(contig_offset[q], x_lo, span)
            lane_bot, _ = lane_y(lane_of[q])
            ax.plot([anchor_x, slot_x], [lane_bot, LBL_TOP],
                    color='#888', linewidth=0.4, zorder=1, clip_on=False)
            name_base = 11 if one_plot else 6
            ax.text(slot_x, LBL_TOP - 0.01, q,
                    ha='left', va='top', rotation=-30,
                    rotation_mode='anchor',
                    fontsize=_name_fontsize(q, name_base) * font_scale,
                    color='#222', clip_on=False)

    # axis cosmetics
    if has_cyt:
        bsp = bands_span(cytobands)
        lbl_lo, lbl_hi = bsp[0], bsp[1]
        lbl_x0, lbl_x1 = rb_x0, rb_x1
    else:
        lbl_lo, lbl_hi = x_lo, x_hi
        lbl_x0, lbl_x1 = 0.0, width_frac
    div, suf = choose_unit([lbl_lo, lbl_hi])
    top_anchor = CYT_TOP if has_cyt else REF_TOP
    ax.text(lbl_x0, top_anchor + 0.01, fmt_unit(lbl_lo, div, suf),
            ha='left',  va='bottom', fontsize=7 * font_scale, color='#555')
    ax.text(lbl_x1, top_anchor + 0.01, fmt_unit(lbl_hi, div, suf),
            ha='right', va='bottom', fontsize=7 * font_scale, color='#555')
    # Empty panels (cytoband-only chroms with no mapped contigs): keep the ref
    # name on the left and draw a dashed open box in the contig lane as a
    # "missing contig" placeholder. No count and no title.
    if not contigs:
        lane_bot, lane_top = _lane_y(MID_TOP, MID_BOT, 1, 0, fill=BAR_FILL)
        ax.add_patch(patches.Rectangle(
            (rb_x0, lane_bot), rb_x1 - rb_x0, lane_top - lane_bot,
            fill=False, edgecolor='#888', linewidth=0.8,
            linestyle=(0, (4, 3)), zorder=3))

    ref_base = 11 if use_left_labels else 8
    ax.text(-0.02, (REF_BOT + REF_TOP) / 2, _fit_name(ref_ctg),
            ha='right', va='center',
            fontsize=_name_fontsize(ref_ctg, ref_base) * font_scale,
            color='#333')
    if contigs:
        count_str = f'{len(contigs)} contigs'
        show_count = len(contigs) > 1
        if use_left_labels and show_count:
            ax.text(-0.02, REF_BOT - 0.005,
                    f'({count_str})',
                    ha='right', va='top',
                    fontsize=9 * font_scale, color='#666')
        if not use_left_labels and not one_plot:
            title = f'{ref_ctg}  ({count_str})' if show_count else ref_ctg
            ax.set_title(title, fontsize=9 * font_scale, pad=15)

    ax.set_xlim(-0.16, 1.06)
    ax.set_ylim(0, 1)
    ax.axis('off')


# ── ctg-centric panel ─────────────────────────────────────────────────────

def draw_ctg_panel(ax, ctg, major_ref, major_blocks, sec_blocks,
                   q_lo, q_hi, global_span=None, font_scale=1.0,
                   centromere=None):
    panel_span = max(q_hi - q_lo, 1)
    span = global_span if global_span else panel_span
    width_frac = panel_span / span

    has_sec = len(sec_blocks) > 0
    if has_sec:
        # Group secondary blocks by t_ctg so blocks of the same non-major ref
        # share a lane and a single label.
        sec_groups = []   # (q_lo_grp, q_hi_grp, t_ctg, sub_df)
        for t_ctg, sub in sec_blocks.groupby('t_ctg', sort=False):
            qlo = int(sub['q_start'].min())
            qhi = int(sub['q_end'].max())
            sec_groups.append((qlo, qhi, t_ctg, sub))
        sec_intervals = [(g[0], g[1], (g[2], g[3])) for g in sec_groups]
        sec_packed, n_sec_lanes = lane_pack(sec_intervals)
    else:
        sec_packed, n_sec_lanes = [], 0

    # bands stacked bottom-up: ctg + gap + maj [+ gap + sec(n_sec) + gap + lbl]
    u, g = float(PANEL_UNIT), float(GAP_UNIT)
    lbl_units = LABEL_BAND_UNITS
    total = u + g + u
    if has_sec:
        total += g + n_sec_lanes * u + g + lbl_units

    CTG_BOT = 0.0
    CTG_TOP = CTG_BOT + u / total
    MID_BOT = CTG_TOP + g / total
    MID_TOP = MID_BOT + u / total
    if has_sec:
        SEC_BOT = MID_TOP + g / total
        SEC_TOP = SEC_BOT + n_sec_lanes * u / total
        LBL_BOT = SEC_TOP + g / total
        LBL_TOP = LBL_BOT + lbl_units / total
    else:
        SEC_BOT = SEC_TOP = MID_TOP
        LBL_BOT = LBL_TOP = MID_TOP

    # The contig bar uses q-coords (normalised so the row width reflects the
    # contig's actual length). The major-ref bar gets its OWN t-coord scale,
    # so each block's drawn width on the mid bar is its true ref length —
    # any q/t length mismatch (insertion/deletion) shows up as a skewed
    # ribbon between the two bars.
    if len(major_blocks):
        t_lo = int(major_blocks['t_start'].min())
        t_hi = int(major_blocks['t_end'].max())
    else:
        t_lo, t_hi = 0, 1
    # extend t-range to include the centromere so it's always visible on the
    # major-ref bar, even when the contig only maps to one arm
    if centromere is not None:
        t_lo = min(t_lo, centromere[0])
        t_hi = max(t_hi, centromere[1])
    t_span = max(t_hi - t_lo, 1)
    mb_width = t_span / span    # same global denominator as width_frac

    # Anchor the major-ref bar on the panel so its smallest-t block lines up
    # horizontally with the corresponding region on the contig bar.
    # mirrors the ref-centric backbone-offset rule.
    if len(major_blocks):
        anchor = major_blocks.loc[major_blocks['t_start'].idxmin()]
        t_anchor = int(anchor['t_start'])
        q_anchor = int(anchor['q_start'])
        ref_x0 = max(0.0,
                     (q_anchor - q_lo) / span - (t_anchor - t_lo) / span)
    else:
        ref_x0 = 0.0

    def nq(x):
        return _norm(x, q_lo, span)

    def nt(x):
        return ref_x0 + (x - t_lo) / span

    # contig bar (thinner than its band)
    cb_bot, cb_top = _shrink(CTG_BOT, CTG_TOP)
    _bar(ax, 0, cb_bot, width_frac, cb_top - cb_bot, CTG_BAR, z=2)
    # major-ref bar (thinner than its band) — width reflects actual t-span
    mb_bot, mb_top = _shrink(MID_BOT, MID_TOP)
    _bar(ax, ref_x0, mb_bot, mb_width, mb_top - mb_bot, REF_BAR, z=2)

    # Centromere on the major-ref bar — t-range is extended above to include
    # it, so this draws unconditionally when the cytoband entry is present.
    if centromere is not None:
        cs, ce = centromere
        if ce > cs:
            x0 = nt(cs)
            x1 = nt(ce)
            mid = (x0 + x1) / 2
            bar_h = mb_top - mb_bot
            tip_y = (mb_bot + mb_top) / 2
            y_top = mb_top + bar_h * 0.25
            y_bot = mb_bot - bar_h * 0.25
            for pts in ([(x0, y_top), (x1, y_top), (mid, tip_y)],
                        [(x0, y_bot), (x1, y_bot), (mid, tip_y)]):
                ax.add_patch(plt.Polygon(pts, closed=True,
                                         facecolor='#333333',
                                         edgecolor='#000000',
                                         linewidth=0.4, zorder=4))

    # major-ref blocks + ribbons between t-coords (mid) and q-coords (contig)
    for _, blk in major_blocks.iterrows():
        color, strand = strand_color(blk['n_pos'], blk['n_neg'])
        if color is None:
            continue
        ts = nt(blk['t_start'])
        te = nt(blk['t_end'])
        qs = nq(blk['q_start'])
        qe = nq(blk['q_end'])
        _bar(ax, ts, mb_bot, te - ts, mb_top - mb_bot,
             color, alpha=0.85, z=3)
        _bar(ax, qs, cb_bot, qe - qs, cb_top - cb_bot,
             color, alpha=0.85, z=3)
        _ribbon(ax,
                ts, te, mb_bot,
                qs, qe, cb_top,
                color, strand)

    # secondary blocks — one lane per non-major t_ctg
    if has_sec:
        for g_lo, g_hi, payload, li in sec_packed:
            t_ctg_g, sub = payload
            y_bot, y_top = _lane_y_up(SEC_BOT, SEC_TOP,
                                      n_sec_lanes, li, fill=BAR_FILL)
            # strike-through across the group's full q-span, matching the
            # ref-centric contig backbone style.
            bb_x0 = _norm(g_lo, q_lo, span)
            bb_x1 = _norm(g_hi, q_lo, span)
            rng = np.random.default_rng(abs(hash(t_ctg_g)) % (2**32))
            lane_h = y_top - y_bot
            jitter = (rng.random() - 0.5) * lane_h * 0.5
            y_strike = (y_bot + y_top) / 2 + jitter
            ax.plot([bb_x0, bb_x1], [y_strike, y_strike],
                    color='black', linewidth=0.6, zorder=10)
            for _, blk in sub.iterrows():
                color, strand = strand_color(blk['n_pos'], blk['n_neg'])
                if color is None:
                    continue
                sx = _norm(int(blk['q_start']), q_lo, span)
                ex = _norm(int(blk['q_end']),   q_lo, span)
                _bar(ax, sx, y_bot, ex - sx, y_top - y_bot,
                     color, alpha=0.85, z=3)
                _bar(ax, sx, cb_bot, ex - sx, cb_top - cb_bot,
                     color, alpha=0.85, z=3)
                _ribbon(ax,
                        sx, ex, y_bot,
                        sx, ex, cb_top,
                        color, strand)

        # secondary-ref name strip above the sec band — one name per t_ctg
        # group, leader line down to the group's leftmost q-position.
        items_sorted = sorted(sec_packed, key=lambda it: it[0])
        n = len(items_sorted)
        for i, (g_lo, g_hi, payload, li) in enumerate(items_sorted):
            t_ctg_g, _sub = payload
            slot_x = (i + 0.5) / n * width_frac if n else 0.0
            sx = _norm(g_lo, q_lo, span)
            _, lane_top = _lane_y_up(SEC_BOT, SEC_TOP,
                                     n_sec_lanes, li, fill=BAR_FILL)
            ax.plot([sx, slot_x], [lane_top, LBL_BOT],
                    color='#888', linewidth=0.4, zorder=1, clip_on=False)
            ax.text(slot_x, LBL_BOT + 0.01, t_ctg_g,
                    ha='left', va='bottom', rotation=30,
                    rotation_mode='anchor',
                    fontsize=5 * font_scale, color='#444', clip_on=False)

    # labels
    div, suf = choose_unit([q_lo, q_hi])
    ax.text(0,          CTG_BOT - 0.02, fmt_unit(q_lo, div, suf),
            ha='left',  va='top', fontsize=7 * font_scale, color='#555')
    ax.text(width_frac, CTG_BOT - 0.02, fmt_unit(q_hi, div, suf),
            ha='right', va='top', fontsize=7 * font_scale, color='#555')
    ax.text(-0.02, (CTG_BOT + CTG_TOP) / 2, _fit_name(ctg),
            ha='right', va='center',
            fontsize=_name_fontsize(ctg, 8) * font_scale, color='#333')
    n_refs = 1 + (sec_blocks['t_ctg'].nunique() if has_sec else 0)
    if n_refs > 1:
        ax.text(-0.02, CTG_BOT - 0.005, f'({n_refs} refs)',
                ha='right', va='top',
                fontsize=7 * font_scale, color='#666')
    ax.text(-0.02, (MID_BOT + MID_TOP) / 2, _fit_name(major_ref),
            ha='right', va='center',
            fontsize=_name_fontsize(major_ref, 8) * font_scale, color='#333')

    ax.set_xlim(-0.14, 1.06)
    ax.set_ylim(0, 1)
    ax.axis('off')


# ── page layout / drivers ─────────────────────────────────────────────────

PANEL_UNIT       = 10  # height of a single lane (arbitrary units)
GAP_UNIT         = 10  # vertical gap between sub-bands inside a subplot
LABEL_BAND_UNITS = 25  # height of the contig-name strip in ref-centric
INTER_SUBPLOT    = 8   # gap between adjacent subplots on a page
ROWS_PER_PAGE    = 5


PAGE_WIDTH_IN  = 16
PAGE_HEIGHT_IN = 9    # fixed 16:9 page ratio


def render_pdf(panels, draw_fn, outfile, title_prefix,
               rows=ROWS_PER_PAGE, one_page=False, one_plot=False,
               show_centromere=False):
    """`one_page`: all panels on a single tall page (height scales with count).
    `one_plot`:   one panel per page (standard 16:9 page each)."""
    ncols = 1
    if one_page:
        ncols = 2
        rows = max(math.ceil(len(panels) / ncols), 1)
    elif one_plot:
        rows = 1
    """One subplot per row on a fixed 16:9 page. All subplots get the same
    height; sub-band proportions inside each subplot still scale with that
    panel's lane counts (handled in the draw functions)."""
    panels_per_page = rows * ncols
    n_pages = max(1, math.ceil(len(panels) / panels_per_page))

    fwd_patch = patches.Patch(color=FWD_COLOR, alpha=0.8, label='Forward')
    rev_patch = patches.Patch(color=REV_COLOR, alpha=0.8, label='Reverse')
    legend_handles = [fwd_patch, rev_patch]

    with PdfPages(outfile) as pdf:
        for page in range(n_pages):
            page_panels = panels[page * panels_per_page:
                                 (page + 1) * panels_per_page]
            hspace = INTER_SUBPLOT / PANEL_UNIT
            if one_page:
                hspace *= 0.3

            page_w = PAGE_WIDTH_IN * ncols
            page_h = (PAGE_HEIGHT_IN * rows / ROWS_PER_PAGE
                      if one_page else PAGE_HEIGHT_IN)
            fig = plt.figure(figsize=(page_w, page_h))
            # Reserve fixed *inches* (not figure fractions) for the title and
            # legend, so very tall one-page figures don't waste enormous
            # space top/bottom.
            top_in, bot_in = 0.85, 0.35
            top_frac = 1.0 - top_in / page_h
            bot_frac = bot_in / page_h
            gs = fig.add_gridspec(nrows=rows, ncols=ncols, hspace=hspace,
                                  wspace=0.12,
                                  top=top_frac, bottom=bot_frac,
                                  left=0.05, right=0.98)
            for slot in range(rows * ncols):
                # column-major fill: col 0 top→bottom, then col 1
                row = slot % rows
                col = slot // rows
                panel_idx = col * rows + row
                ax = fig.add_subplot(gs[row, col])
                if panel_idx >= len(page_panels):
                    ax.axis('off')
                    continue
                draw_fn(ax, page_panels[panel_idx])

            fig.legend(handles=legend_handles,
                       loc='lower center', ncol=len(legend_handles),
                       fontsize=11, frameon=False,
                       bbox_to_anchor=(0.5, bot_frac * 0.25))
            fig.suptitle(f'{title_prefix}  (page {page + 1}/{n_pages})',
                         fontsize=13, y=1.0 - 0.2 / page_h)
            pdf.savefig(fig)
            plt.close(fig)
            print(f'  page {page + 1}/{n_pages} written')


def run_ref_centric(df, major, outfile, one_page=False, one_plot=False,
                    cytobands=None):
    """One panel per reference. Each panel includes EVERY contig that has at
    least one block on that reference (not just contigs whose major ref it
    is), so a contig mapping to multiple refs appears on multiple panels.
    `major` is still used to mark each placement as primary vs shared."""
    contigs_by_ref = {}
    for q, t in zip(df['query_ctg'], df['t_ctg']):
        contigs_by_ref.setdefault(t, set()).add(q)

    # Add empty entries for cytoband chroms with no mapped contigs, so the
    # ideogram still gets a panel. parse_cytobands indexes the same list under
    # both the full chrom name and its 'chrN' alias — dedup by list identity
    # and prefer the longer (more specific) key.
    if cytobands:
        seen_ids = set()
        unique_keys = {}
        for k, v in cytobands.items():
            i = id(v)
            if i in seen_ids:
                if len(k) > len(unique_keys[i]):
                    unique_keys[i] = k
                continue
            seen_ids.add(i)
            unique_keys[i] = k
        for k in unique_keys.values():
            contigs_by_ref.setdefault(k, set())

    panels = []
    for ref_ctg in sorted(contigs_by_ref.keys(), key=natural_key):
        ctgs = contigs_by_ref[ref_ctg]
        # order contigs by their start on this ref (stable)
        ctgs_sorted = sorted(
            ctgs,
            key=lambda q: int(df[(df['query_ctg'] == q) &
                                 (df['t_ctg'] == ref_ctg)]['t_start'].min())
        ) if ctgs else []
        ref_blocks_by_ctg = {}
        sec_blocks_by_ctg = {}
        primary_of = {}
        x_lo = math.inf
        x_hi = -math.inf
        for q in ctgs_sorted:
            qb = df[df['query_ctg'] == q]
            mb = qb[qb['t_ctg'] == ref_ctg]
            sb = qb[qb['t_ctg'] != ref_ctg]
            ref_blocks_by_ctg[q] = mb
            sec_blocks_by_ctg[q] = sb
            primary_of[q] = (major.get(q) == ref_ctg)
            offset = int(mb['t_start'].min())
            x_lo = min(x_lo, offset)
            x_hi = max(x_hi, int(mb['t_end'].max()))
            # contig-lane extent: backbone runs offset → offset + ctg_length
            q_lo = int(mb['q_start'].min())
            q_hi = int(mb['q_end'].max())
            if len(sb):
                q_lo = min(q_lo, int(sb['q_start'].min()))
                q_hi = max(q_hi, int(sb['q_end'].max()))
            x_hi = max(x_hi, offset + (q_hi - q_lo))

        # estimate contig-lane count to size the subplot row height
        intervals = []
        for q in ctgs_sorted:
            mb = ref_blocks_by_ctg[q]
            sb = sec_blocks_by_ctg[q]
            off = int(mb['t_start'].min())
            qmin = int(mb['q_start'].min())
            qmax = int(mb['q_end'].max())
            if len(sb):
                qmin = min(qmin, int(sb['q_start'].min()))
                qmax = max(qmax, int(sb['q_end'].max()))
            intervals.append((off, off + (qmax - qmin), q))
        _, n_mid_lanes = lane_pack(intervals)
        total_units = PANEL_UNIT + GAP_UNIT + n_mid_lanes * PANEL_UNIT

        # top-kmer contig on this ref (used as the lone label under --one-page);
        # ranked by on-panel kmer only so a "shared" contig with little signal
        # here doesn't outrank a primary one.
        top_ctg = max(
            ctgs_sorted,
            key=lambda q: int(ref_blocks_by_ctg[q]['kmer'].sum())
        ) if ctgs_sorted else None

        bands = lookup_cytobands(cytobands, ref_ctg) if cytobands else None
        cen = centromere_span(bands)
        # If cytobands are present, extend to the full chromosome so the
        # ideogram is shown end-to-end (same spirit as the prior centromere
        # extension, just broader).
        bspan = bands_span(bands)
        if bspan is not None:
            x_lo = min(x_lo, bspan[0])
            x_hi = max(x_hi, bspan[1])
        elif cen is not None:
            x_lo = min(x_lo, cen[0])
            x_hi = max(x_hi, cen[1])

        panels.append({
            'ref_ctg': ref_ctg,
            'ref_blocks_by_ctg': ref_blocks_by_ctg,
            'sec_blocks_by_ctg': sec_blocks_by_ctg,
            'primary_of': primary_of,
            'x_lo': x_lo, 'x_hi': x_hi,
            'total_units': total_units,
            'top_ctg': top_ctg,
            'centromere': cen,
            'cytobands': bands,
        })

    n_prim = sum(v for p in panels for v in p['primary_of'].values())
    n_shar = sum(1 - v for p in panels for v in p['primary_of'].values())
    print(f'ref-centric: {len(panels)} reference panels '
          f'({n_prim} primary + {n_shar} shared contig placements)')
    global_span = max((p['x_hi'] - p['x_lo'] for p in panels), default=1)

    fscale = 1.5 if one_plot else 1.0

    def _draw(ax, p):
        labeled = [p['top_ctg']] if one_page else None
        draw_ref_panel(ax, p['ref_ctg'],
                       p['ref_blocks_by_ctg'], p['sec_blocks_by_ctg'],
                       p['x_lo'], p['x_hi'],
                       global_span=global_span,
                       labeled_contigs=labeled,
                       font_scale=fscale,
                       one_plot=one_plot,
                       centromere=p.get('centromere'),
                       primary_of=p['primary_of'],
                       cytobands=p.get('cytobands'))

    render_pdf(panels, _draw, outfile,
               title_prefix='Reference-centric block mapping',
               one_page=one_page, one_plot=one_plot,
               show_centromere=bool(cytobands))


def run_ctg_centric(df, major, outfile, one_page=False, one_plot=False,
                    cytobands=None):
    if one_page:
        print('  hint: --one-page with ctg-centric can be very dense; '
              'consider --min-block-kmer / --min-pair-kmer to thin labels')
    panels = []
    contigs = sorted(major.keys(), key=natural_key)
    for q in contigs:
        major_ref = major[q]
        qb = df[df['query_ctg'] == q]
        mb = qb[qb['t_ctg'] == major_ref]
        sb = qb[qb['t_ctg'] != major_ref]
        q_lo = int(qb['q_start'].min())
        q_hi = int(qb['q_end'].max())
        if len(sb):
            sec_intervals = [(int(sub['q_start'].min()),
                              int(sub['q_end'].max()), None)
                             for _, sub in sb.groupby('t_ctg', sort=False)]
            n_sec_lanes = lane_pack(sec_intervals)[1]
        else:
            n_sec_lanes = 0
        total_units = PANEL_UNIT + GAP_UNIT + PANEL_UNIT
        if n_sec_lanes > 0:
            total_units += GAP_UNIT + n_sec_lanes * PANEL_UNIT

        cen = (centromere_span(lookup_cytobands(cytobands, major_ref))
               if cytobands else None)

        panels.append({
            'ctg': q, 'major_ref': major_ref,
            'major_blocks': mb, 'sec_blocks': sb,
            'q_lo': q_lo, 'q_hi': q_hi,
            'total_units': total_units,
            'centromere': cen,
        })

    print(f'ctg-centric: {len(panels)} contig panels')
    global_span = max((p['q_hi'] - p['q_lo'] for p in panels), default=1)

    fscale = 1.5 if one_plot else 1.0

    def _draw(ax, p):
        draw_ctg_panel(ax, p['ctg'], p['major_ref'],
                       p['major_blocks'], p['sec_blocks'],
                       p['q_lo'], p['q_hi'],
                       global_span=global_span,
                       font_scale=fscale,
                       centromere=p.get('centromere'))

    render_pdf(panels, _draw, outfile,
               title_prefix='Contig-centric block mapping',
               one_page=one_page, one_plot=one_plot,
               show_centromere=bool(cytobands))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input')
    ap.add_argument('output', nargs='?', default='m_panel_plot.pdf')
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--ref-centric', action='store_true')
    mode.add_argument('--ctg-centric', action='store_true')
    ap.add_argument('--min-block-kmer', type=int, default=0,
                    help='drop blocks with (n_pos+n_neg) below this (default 0)')
    ap.add_argument('--min-pair-kmer', type=int, default=0,
                    help='drop a (contig, ref) pair if its SUMMED (n_pos+n_neg) '
                         'across all blocks is below this (default 0)')
    ap.add_argument('--strand', default=None,
                    help='per-(contig, ref) strand file (col1=ctg, col2=ref, '
                         'col3=tau). When the (ctg, major_ref) tau is < 0, '
                         'the contig is reverse-complemented: q_new = L - q_old '
                         '(L = max q_end of the contig) and n_pos/n_neg are '
                         'swapped. All of that contig\'s blocks flip together. '
                         'Pairs missing from the file are left as-is.')
    ap.add_argument('--cytobands', default=None,
                    help='UCSC cytoband BED (chrom start end name stain). '
                         'In --ref-centric, each panel gets a full ideogram '
                         'strip below the reference bar (panel x-range is '
                         'extended end-to-end across the chrom). The '
                         'centromere (acen bands) is also overlaid as an '
                         'hourglass on the ref bar in --ref-centric and on '
                         'the major-ref bar in --ctg-centric.')
    layout = ap.add_mutually_exclusive_group()
    layout.add_argument('--one-page', action='store_true',
                        help='put every panel on a single tall page; '
                             'for ref-centric, one contig per lane is labeled '
                             '(the highest-kmer contig in each lane)')
    layout.add_argument('--one-plot', action='store_true',
                        help='one panel per page (full detail, all labels)')
    args = ap.parse_args()

    df = parse_blocks(args.input)
    print(f'parsed {len(df)} block rows')
    df, major = assign_major_ref(df,
                                 min_kmer=args.min_block_kmer,
                                 min_pair_kmer=args.min_pair_kmer)
    print(f'after filter: {len(df)} rows; {len(major)} contigs assigned a major ref')

    if args.strand:
        strand = parse_strand_tsv(args.strand)
        print(f'parsed {len(strand)} (ctg, ref) strand rows from '
              f'{args.strand}')
        df, flipped = apply_strand_flips(df, major, strand)
        print(f'reversed {len(flipped)} contigs based on (ctg, major_ref) '
              f'col10 > col9')

    cytobands = parse_cytobands(args.cytobands) if args.cytobands else None
    if cytobands:
        n_bands = sum(len(v) for v in cytobands.values())
        print(f'parsed {n_bands} cytoband rows across {len(cytobands)} '
              f'chrom keys from {args.cytobands}')

    if args.ref_centric:
        run_ref_centric(df, major, args.output,
                        one_page=args.one_page, one_plot=args.one_plot,
                        cytobands=cytobands)
    else:
        run_ctg_centric(df, major, args.output,
                        one_page=args.one_page, one_plot=args.one_plot,
                        cytobands=cytobands)

    print(f'\ndone. output: {args.output}')


if __name__ == '__main__':
    main()
