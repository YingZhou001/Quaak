# quick_block

k-mer based synteny block detection and structural variant analysis between two genome assemblies.

## Overview

`quick_block` maps unique k-mers from a reference assembly onto a target assembly, groups collinear k-mer matches into synteny blocks, and identifies regions of the reference that are absent in the target.

## Build

Requires zlib and a C99 compiler.

```
gcc -O2 -o quick_block quick_block.c -lz -lm
```

`khashl.h` must be present in the same directory.

## Usage

```
quick_block [--cut gap_skip_cut,gap_diff_cut] <ref.gz> <tar.gz> <outpref>
```

| Argument | Description |
|---|---|
| `ref.gz` | k-mer BED file for the reference assembly (gzip-compressed) |
| `tar.gz` | k-mer BED file for the target assembly (gzip-compressed) |
| `outpref` | Prefix for all output files |

### Options

| Option | Default | Description |
|---|---|---|
| `--cut skip,diff` | `100000000,100000000` | Gap thresholds for block splitting. `skip` is the maximum allowed gap between consecutive k-mers on either sequence; `diff` is the maximum allowed difference between the two gap sizes. |

## Input format

Both `ref.gz` and `tar.gz` are tab/space-separated files with one k-mer per line:

```
<contig>  <start>  <end>  <kmer>  <misc>  <strand>
```

- Coordinates are 0-based.
- `<kmer>` may contain a `:` suffix (e.g. `ACGT:count`) which is stripped before use.
- `<strand>` is `+` or `-`.
- K-mers appearing more than once in the reference are discarded.

## Output files

### `<outpref>.ctg-strand.csv`

Per-contig synteny summary. One row per (target contig, reference contig) pair. Tab-separated, no header.

| Column | Description |
|---|---|
| `ctg` | Target contig name |
| `tar_fro` | Start of matched region on target |
| `tar_to` | End of matched region on target |
| `ref_ctg` | Reference contig name |
| `ref_fro` | Start of matched region on reference |
| `ref_to` | End of matched region on reference |
| `tau` | Kendall's tau rank correlation of k-mer positions |
| `pvalue` | Two-tailed p-value for tau (normal approximation) |
| `n_plus` | Number of k-mers on `+` strand |
| `n_minus` | Number of k-mers on `-` strand |

### `<outpref>.block.csv.gz`

Collinear synteny blocks (gzip-compressed). One row per block. Tab-separated, no header.

| Column | Description |
|---|---|
| `ctg` | Target contig name |
| `block_idx` | Block index (resets to 0 per target contig) |
| `tar_fro` | Block start on target |
| `tar_to` | Block end on target |
| `ref_ctg` | Reference contig name |
| `ref_fro` | Block start on reference |
| `ref_to` | Block end on reference |
| `n_plus` | K-mers on `+` strand |
| `n_minus` | K-mers on `-` strand |
| `kmers` | First and last k-mer of the block, comma-separated |

### `<outpref>.kmer-summary.csv`

CNV pivot table. Rows are reference contigs (in input order); columns are copy-number values (0, 1, 2, …). Each cell is the number of k-mers from that contig with that copy number in the target. Tab-separated with a header row.

### `<outpref>.missing.csv`

Regions of the reference with no k-mer support in the target (copy number = 0). One row per contiguous run of zero-count k-mers. Tab-separated, no header.

| Column | Description |
|---|---|
| `ctg` | Reference contig name |
| `fro` | Start of missing region |
| `to` | End of missing region |
| `in_kmers` | First and last k-mer inside the missing region, comma-separated |
| `left_pos` | End coordinate of the flanking block to the left (`NA` if none) |
| `right_pos` | End coordinate of the flanking block to the right (`NA` if none) |
| `out_kmers` | Last k-mer of left flank and last k-mer of right flank, comma-separated |
| `n` | Number of k-mers in the missing region |

## Algorithm

1. **Load reference** — parse `ref.gz`, build a hash table of unique k-mers keyed by sequence. Duplicate k-mers (appearing more than once) are removed.
2. **Process target** — parse `tar.gz`, look up each k-mer in the reference table. Matched k-mers are buffered per target contig. When the contig changes, `block_split` is called.
3. **Block splitting** — within a target contig buffer, consecutive k-mer matches are merged into a block if they share the same reference contig, strand, order, and both inter-k-mer gaps are within the cut thresholds. Kendall's tau is computed for each (target contig, reference contig) pair using an O(n log n) merge-sort inversion count.
4. **CNV summary** — after target processing, each reference k-mer carries a hit count. These counts are aggregated into the pivot table.
5. **Missing region detection** — reference k-mers are sorted by position within each contig and run-length encoded by count. Runs with count = 0 are written to `missing.csv`.
