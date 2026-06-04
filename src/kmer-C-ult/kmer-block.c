/*
 * kmer-block.c  –  C port of kmer-block.py
 *
 * Usage: kmer-block <ref.gz> <tar.gz> <outpref>
 *
 * Outputs:
 *   <outpref>.ctg-strand.tsv     – per-contig synteny summary (Kendall tau)
 *   <outpref>.block.tsv.gz       – synteny blocks
 *   <outpref>.kmer-summary.tsv   – k-mer CNV pivot table
 *   <outpref>.cnv.tsv.gz         – all CNV (run-length-encoded count) blocks
 *
 * Compile:
 *   gcc -O2 -o kmer-block kmer-block.c -lz -lm
 */

#define _POSIX_C_SOURCE 200809L   /* expose strdup, strnlen under -std=c99 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <zlib.h>
#include "khashl.h"

#define MAX_STR  256
#define MAX_PATH 512
#define LINE_BUF 2048

/* ------------------------------------------------------------------ */
/* Hash table types                                                    */
/* ------------------------------------------------------------------ */

/* string key -> void* value  (generic pointer map, cast at use site) */
KHASHL_CMAP_INIT(KH_LOCAL, kh_strmap_t, kh_strmap, const char*, void*,
                 kh_hash_str, kh_eq_str)

/* string key -> int value  (used for kmer counts and ctg counts)     */
KHASHL_CMAP_INIT(KH_LOCAL, kh_strint_t, kh_strint, const char*, int,
                 kh_hash_str, kh_eq_str)

/* uint key -> void* value  (used for CNV cp -> inner map)            */
KHASHL_MAP_INIT(KH_LOCAL, kh_intmap_t, kh_intmap, khint_t, void*,
                kh_hash_uint32, kh_eq_generic)

/* ------------------------------------------------------------------ */
/* Data structures                                                     */
/* ------------------------------------------------------------------ */

typedef struct {
    char kmer[MAX_STR];   /* also used as hash key (pointer to this field) */
    char ctg [MAX_STR];
    long fro, to;
    char strand;
    int  count;           /* -1 = duplicate; >=0 = target hit count        */
} KmerEntry;

typedef struct {
    char ctg    [MAX_STR];
    long fro, to;
    char ref_ctg[MAX_STR];
    long ref_fro, ref_to;
    char kmer   [MAX_STR];
    char strand;
} DARow;

typedef struct {
    char  ref_ctg[MAX_STR];   /* used as hash key */
    long *tar_pos, *ref_pos;
    char *strands;
    int   size, cap;
    long  tar_fro, tar_to, ref_fro, ref_to;
} RefCtgAcc;

typedef struct {
    char   ctg    [MAX_STR];
    long   tar_fro, tar_to;
    char   ref_ctg[MAX_STR];
    long   ref_fro, ref_to;
    double coef, pvalue;
    int    n_plus, n_minus;
} SummaryRow;

typedef struct {
    DARow *rows;
    int    size, cap;
} Block;

typedef struct {
    Block *data;
    int    size, cap;
} BlockArr;

typedef struct {
    long fro, to;
    int  count;
    char kmer[MAX_STR];
} KmerPos;

typedef struct {
    char    ctg [MAX_STR];   /* used as hash key */
    KmerPos *data;
    int      size, cap;
} CtgKmerList;

typedef struct {
    char ctg       [MAX_STR];
    char kmer_first[MAX_STR];
    char kmer_last [MAX_STR];
    long fro, to;
    int  count, n;
} CnvBlock;

/* ------------------------------------------------------------------ */
/* Global reference k-mer table                                        */
/* ------------------------------------------------------------------ */
static kh_strmap_t *ref_dict = NULL;

/* ------------------------------------------------------------------ */
/* String helpers                                                      */
/* ------------------------------------------------------------------ */
/* Bounded string copy that always NUL-terminates within MAX_STR bytes.
 * Plain strncpy(dst, src, MAX_STR-1) does *not* terminate when src is
 * 255 chars long (sscanf %255s produces exactly that), leaving the
 * trailing byte whatever malloc handed back. */
static void copy_str(char *dst, const char *src)
{
    size_t n = strnlen(src, MAX_STR - 1);
    memcpy(dst, src, n);
    dst[n] = '\0';
}

/* ------------------------------------------------------------------ */
/* Small dynamic-array helpers                                         */
/* ------------------------------------------------------------------ */
static void acc_push(RefCtgAcc *a, long tp, long rp, char s)
{
    if (a->size == a->cap) {
        a->cap   *= 2;
        a->tar_pos = realloc(a->tar_pos, a->cap * sizeof(long));
        a->ref_pos = realloc(a->ref_pos, a->cap * sizeof(long));
        a->strands = realloc(a->strands, a->cap * sizeof(char));
    }
    a->tar_pos[a->size] = tp;
    a->ref_pos[a->size] = rp;
    a->strands[a->size] = s;
    a->size++;
}

static void block_push(Block *b, const DARow *r)
{
    if (b->size == b->cap) {
        b->cap *= 2;
        b->rows = realloc(b->rows, b->cap * sizeof(DARow));
    }
    b->rows[b->size++] = *r;
}

static void ba_push(BlockArr *ba, const Block *b)
{
    if (ba->size == ba->cap) {
        ba->cap *= 2;
        ba->data = realloc(ba->data, ba->cap * sizeof(Block));
    }
    ba->data[ba->size++] = *b;
}

static void free_blocks(BlockArr *ba)
{
    for (int i = 0; i < ba->size; i++) free(ba->data[i].rows);
    free(ba->data);
}

/* ------------------------------------------------------------------ */
/* Kendall's tau  (O(n log n) via merge-sort inversion count)         */
/* ------------------------------------------------------------------ */
typedef struct { long x, y; } XYPair;

static int cmp_xy(const void *a, const void *b)
{
    long dx = ((const XYPair *)a)->x - ((const XYPair *)b)->x;
    if (dx != 0) return dx < 0 ? -1 : 1;
    long dy = ((const XYPair *)a)->y - ((const XYPair *)b)->y;
    return dy < 0 ? -1 : dy > 0 ? 1 : 0;
}

/* Merge-sort arr[0..n-1] in-place using tmp[] as workspace.
 * Returns the number of inversions (= discordant pairs).
 * Equal elements are taken from the left half first, so ties in y
 * are not counted as inversions. */
static long merge_inv(long *arr, long *tmp, int n)
{
    if (n <= 1) return 0;
    int  mid = n / 2;
    long inv = merge_inv(arr, tmp, mid) + merge_inv(arr + mid, tmp + mid, n - mid);
    int i = 0, j = mid, k = 0;
    while (i < mid && j < n) {
        if (arr[i] <= arr[j]) tmp[k++] = arr[i++];
        else { inv += mid - i; tmp[k++] = arr[j++]; }
    }
    while (i < mid) tmp[k++] = arr[i++];
    while (j < n)   tmp[k++] = arr[j++];
    memcpy(arr, tmp, n * sizeof(long));
    return inv;
}

static void cal_tau(const long *x, const long *y, int n,
                    double *tau_out, double *pval_out)
{
    if (n < 3) { *tau_out = NAN; *pval_out = NAN; return; }

    /* sort pairs by x (ties broken by y) so inversions in y == discordant pairs */
    XYPair *pairs = malloc(n * sizeof(XYPair));
    for (int i = 0; i < n; i++) { pairs[i].x = x[i]; pairs[i].y = y[i]; }
    qsort(pairs, n, sizeof(XYPair), cmp_xy);

    long *yarr = malloc(n * sizeof(long));
    long *tmp  = malloc(n * sizeof(long));
    for (int i = 0; i < n; i++) yarr[i] = pairs[i].y;
    free(pairs);

    long   total = (long)n * (n - 1) / 2;
    long   disc  = merge_inv(yarr, tmp, n);
    long   conc  = total - disc;
    double tau   = (double)(conc - disc) / (double)total;
    /* normal approximation: var(tau) = 2(2n+5) / (9n(n-1)).
     * Force double arithmetic so 2*(2n+5) doesn't overflow int for large n. */
    double var   = 2.0*(2.0*n + 5.0) / (9.0*n*(n-1));
    double z     = tau / sqrt(var);
    *tau_out  = tau;
    *pval_out = erfc(fabs(z) / sqrt(2.0));

    free(yarr);
    free(tmp);
}

/* ------------------------------------------------------------------ */
/* block_split                                                         */
/*   DA[0..n-1] are all for the same target contig (sorted by fro).   */
/*   Returns summary rows (one per ref-ctg) and synteny blocks.       */
/* ------------------------------------------------------------------ */
static void block_split(const DARow *DA, int n,
                        long gap_skip_cut, long gap_diff_cut,
                        SummaryRow **summary_out, int *summary_n_out,
                        BlockArr   *blocks_out)
{
    khint_t k;
    int     i, absent;

    /* --- count each kmer within this buffer --- */
    kh_strint_t *kmer_cnt = kh_strint_init();
    for (i = 0; i < n; i++) {
        k = kh_strint_put(kmer_cnt, DA[i].kmer, &absent);
        if (absent) kh_val(kmer_cnt, k) = 0;
        kh_val(kmer_cnt, k)++;
    }

    /* --- accumulate per-ref-ctg position lists and bounding boxes --- */
    kh_strmap_t *acc_map = kh_strmap_init();
    for (i = 0; i < n; i++) {
        const DARow *r = &DA[i];
        k = kh_strmap_get(acc_map, r->ref_ctg);
        RefCtgAcc *acc;
        if (k == kh_end(acc_map)) {
            acc = calloc(1, sizeof(RefCtgAcc));
            copy_str(acc->ref_ctg, r->ref_ctg);
            acc->cap     = 16;
            acc->tar_pos = malloc(acc->cap * sizeof(long));
            acc->ref_pos = malloc(acc->cap * sizeof(long));
            acc->strands = malloc(acc->cap);
            acc->tar_fro = r->fro;     acc->tar_to  = r->to;
            acc->ref_fro = r->ref_fro; acc->ref_to  = r->ref_to;
            k = kh_strmap_put(acc_map, acc->ref_ctg, &absent);
            kh_val(acc_map, k) = acc;
        } else {
            acc = (RefCtgAcc *)kh_val(acc_map, k);
        }
        acc_push(acc, r->fro, r->ref_fro, r->strand);
        if (r->fro     < acc->tar_fro) acc->tar_fro = r->fro;
        if (r->to      > acc->tar_to)  acc->tar_to  = r->to;
        if (r->ref_fro < acc->ref_fro) acc->ref_fro = r->ref_fro;
        if (r->ref_to  > acc->ref_to)  acc->ref_to  = r->ref_to;
    }

    /* --- second pass: group rows into collinear blocks --- */
    blocks_out->data = malloc(16 * sizeof(Block));
    blocks_out->size = 0;
    blocks_out->cap  = 16;
    Block cur = { malloc(16 * sizeof(DARow)), 0, 16 };

    for (i = 0; i < n; i++) {
        const DARow *row = &DA[i];
        if (cur.size == 0) { block_push(&cur, row); continue; }

        const DARow *prev = &cur.rows[cur.size - 1];

        khint_t kn = kh_strint_get(kmer_cnt, row->kmer);
        khint_t ko = kh_strint_get(kmer_cnt, prev->kmer);
        int is_same_kmer_cnt = (kh_val(kmer_cnt, kn) == kh_val(kmer_cnt, ko));
        int is_same_chr      = (strcmp(row->ref_ctg, prev->ref_ctg) == 0);
        int is_same_strand   = (row->strand == prev->strand);

        long gap_tar   = (long)row->fro     - prev->fro;
        long gap_ref   = (long)row->ref_fro - prev->ref_fro;

        if (gap_tar <= 0) {
            fprintf(stderr,
                    "Error: duplicated or unsorted rows in block_split "
                    "(ctg=%s fro=%ld)\n", row->ctg, row->fro);
            exit(1);
        }

        int is_same_order;
        if      (gap_ref == 0)       is_same_order = 0;
        else if (row->strand == '+') is_same_order = (gap_ref > 0);
        else                         is_same_order = (gap_ref < 0);

        long gap_ref_l  = labs(gap_ref);
        int  is_small   = (gap_tar < gap_skip_cut) && (gap_ref_l < gap_skip_cut);
        int  is_sim_gap = labs(gap_ref_l - gap_tar) < gap_diff_cut;

        if (is_same_chr && is_same_strand && is_same_order &&
            is_small && is_sim_gap && is_same_kmer_cnt) {
            block_push(&cur, row);
        } else {
            ba_push(blocks_out, &cur);
            cur = (Block){ malloc(16 * sizeof(DARow)), 0, 16 };
            block_push(&cur, row);
        }
    }
    if (cur.size > 0) ba_push(blocks_out, &cur);
    else              free(cur.rows);

    /* --- build summary (one entry per ref-ctg) --- */
    int nrefs = (int)kh_size(acc_map);
    SummaryRow *summary = malloc(nrefs * sizeof(SummaryRow));
    int si = 0;
    const char *tar_ctg = DA[n-1].ctg;   /* all rows share the same ctg */

    for (khint_t ki = 0; ki < kh_end(acc_map); ki++) {
        if (!kh_exist(acc_map, ki)) continue;
        RefCtgAcc *acc = (RefCtgAcc *)kh_val(acc_map, ki);
        double tau, pval;
        cal_tau(acc->tar_pos, acc->ref_pos, acc->size, &tau, &pval);
        int n_plus = 0;
        for (int j = 0; j < acc->size; j++)
            if (acc->strands[j] == '+') n_plus++;
        SummaryRow *s = &summary[si++];
        copy_str(s->ctg,     tar_ctg);
        s->tar_fro = acc->tar_fro; s->tar_to  = acc->tar_to;
        copy_str(s->ref_ctg, acc->ref_ctg);
        s->ref_fro = acc->ref_fro; s->ref_to  = acc->ref_to;
        s->coef    = tau;           s->pvalue  = pval;
        s->n_plus  = n_plus;        s->n_minus = acc->size - n_plus;
        free(acc->tar_pos); free(acc->ref_pos); free(acc->strands); free(acc);
    }
    kh_strmap_destroy(acc_map);
    *summary_out   = summary;
    *summary_n_out = si;

    kh_strint_destroy(kmer_cnt);
}

/* ------------------------------------------------------------------ */
/* Output helpers                                                      */
/* ------------------------------------------------------------------ */
static void print_summary_line(FILE *fp, const SummaryRow *s, int n)
{
    for (int i = 0; i < n; i++)
        fprintf(fp, "%s\t%ld\t%ld\t%s\t%ld\t%ld\t%.6f\t%.6e\t%d\t%d\n",
                s[i].ctg, s[i].tar_fro, s[i].tar_to,
                s[i].ref_ctg, s[i].ref_fro, s[i].ref_to,
                s[i].coef, s[i].pvalue, s[i].n_plus, s[i].n_minus);
}

/* block_idx resets to 0 per call, matching Python's local variable */
static void print_blocks(gzFile fp, const BlockArr *ba)
{
    for (int bi = 0; bi < ba->size; bi++) {
        const Block *blk  = &ba->data[bi];
        int          nb   = blk->size;
        const DARow *row0 = &blk->rows[0];
        const DARow *row1 = &blk->rows[nb - 1];
        long tar_fro = row0->fro     < row1->fro     ? row0->fro     : row1->fro;
        long tar_to  = row0->to      > row1->to      ? row0->to      : row1->to;
        long ref_fro = row0->ref_fro < row1->ref_fro ? row0->ref_fro : row1->ref_fro;
        long ref_to  = row0->ref_to  > row1->ref_to  ? row0->ref_to  : row1->ref_to;
        int n_plus  = 0;
        for (int j = 0; j < nb; j++)
            if (blk->rows[j].strand == '+') n_plus++;
        int  n_minus = nb - n_plus;
        char kmers[MAX_STR * 2 + 2];
        snprintf(kmers, sizeof(kmers), "%s,%s", row0->kmer, row1->kmer);
        gzprintf(fp, "%s\t%d\t%ld\t%ld\t%s\t%ld\t%ld\t%d\t%d\t%s\n",
                 row0->ctg, bi, tar_fro, tar_to,
                 row0->ref_ctg, ref_fro, ref_to,
                 n_plus, n_minus, kmers);
    }
}

/* ------------------------------------------------------------------ */
/* Comparators                                                         */
/* ------------------------------------------------------------------ */
static int cmp_kmer_pos(const void *a, const void *b)
{
    long d = ((const KmerPos *)a)->fro - ((const KmerPos *)b)->fro;
    return d < 0 ? -1 : d > 0 ? 1 : 0;
}

static int cmp_int(const void *a, const void *b)
{
    khint_t x = *(const khint_t *)a, y = *(const khint_t *)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

static int cmp_strptr(const void *a, const void *b)
{
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}


/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */
int main(int argc, char *argv[])
{
    setvbuf(stdout, NULL, _IOLBF, 0);
    setvbuf(stderr, NULL, _IOLBF, 0);

    long gap_skip_cut = 100000000L;
    long gap_diff_cut = 100000000L;

    /* parse optional --cut skip,diff before positional args */
    int argi = 1;
    while (argi < argc && argv[argi][0] == '-') {
        if (strcmp(argv[argi], "--cut") == 0) {
            if (argi + 1 >= argc) {
                fprintf(stderr, "Error: --cut requires an argument\n");
                return 1;
            }
            if (sscanf(argv[++argi], "%ld,%ld", &gap_skip_cut, &gap_diff_cut) != 2) {
                fprintf(stderr, "Error: --cut expects two comma-separated integers (e.g. 100000,100000)\n");
                return 1;
            }
        } else {
            fprintf(stderr, "Error: unknown option '%s'\n", argv[argi]);
            return 1;
        }
        argi++;
    }

    if (argc - argi != 3) {
        fprintf(stderr, "Usage: kmer-block [--cut gap_skip_cut,gap_diff_cut] <ref.path.gz> <tar.path.gz> <outpref>\n");
        return 1;
    }

    const char *reffile = argv[argi];
    const char *tarfile = argv[argi + 1];
    const char *outpref = argv[argi + 2];

    char out_ctg  [MAX_PATH], out_kmer [MAX_PATH];
    char out_del  [MAX_PATH], out_block[MAX_PATH];
    snprintf(out_ctg,   sizeof(out_ctg),   "%s.ctg-strand.tsv",   outpref);
    snprintf(out_kmer,  sizeof(out_kmer),  "%s.kmer-summary.tsv", outpref);
    snprintf(out_del,   sizeof(out_del),   "%s.cnv.tsv.gz",       outpref);
    snprintf(out_block, sizeof(out_block), "%s.block.tsv.gz",     outpref);

    /* ============================================================
     * 1. Load reference k-mers
     * ============================================================ */
    fprintf(stderr, "[msg] load reference from %s\n", reffile);
    gzFile fp = gzopen(reffile, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", reffile); return 1; }

    ref_dict = kh_strmap_init();
    char line[LINE_BUF];

    /* track ref-ctg first-appearance order */
    kh_strint_t *ctg_seen    = kh_strint_init();
    int          n_ctg_ord   = 0, cap_ctg_ord = 64;
    char       **ctg_ord     = malloc(cap_ctg_ord * sizeof(char *));

    while (gzgets(fp, line, sizeof(line))) {
        char ctg[MAX_STR], kmer[MAX_STR], dummy[MAX_STR], strand_c;
        long fro, to;
        if (sscanf(line, "%255s %ld %ld %255s %255s %c",
                   ctg, &fro, &to, kmer, dummy, &strand_c) != 6) continue;
        char *colon = strchr(kmer, ':');
        if (colon) *colon = '\0';

        khint_t k = kh_strmap_get(ref_dict, kmer);
        if (k == kh_end(ref_dict)) {
            KmerEntry *e = calloc(1, sizeof(KmerEntry));
            copy_str(e->kmer,   kmer);
            copy_str(e->ctg,    ctg);
            e->fro    = fro; e->to = to;
            e->strand = strand_c;
            e->count  = 0;
            int absent;
            k = kh_strmap_put(ref_dict, e->kmer, &absent);
            kh_val(ref_dict, k) = e;

            /* record ctg on first appearance among newly inserted kmers */
            int absent2;
            khint_t kc = kh_strint_put(ctg_seen, ctg, &absent2);
            if (absent2) {
                char *sc = strdup(ctg);
                kh_key(ctg_seen, kc) = sc;
                if (n_ctg_ord == cap_ctg_ord) {
                    cap_ctg_ord *= 2;
                    ctg_ord = realloc(ctg_ord, cap_ctg_ord * sizeof(char *));
                }
                ctg_ord[n_ctg_ord++] = sc;
            }
        } else {
            ((KmerEntry *)kh_val(ref_dict, k))->count = -1; /* duplicate */
        }
    }
    gzclose(fp);
    kh_strint_destroy(ctg_seen);

    /* remove duplicates */
    int removed = 0;
    for (khint_t ki = 0; ki < kh_end(ref_dict); ki++) {
        if (!kh_exist(ref_dict, ki)) continue;
        KmerEntry *e = (KmerEntry *)kh_val(ref_dict, ki);
        if (e->count == -1) {
            kh_strmap_del(ref_dict, ki);
            free(e);
            removed++;
        }
    }
    fprintf(stderr, "[msg] %d kmer loaded, %d removed due to duplication in ref\n",
            (int)kh_size(ref_dict), removed);

    /* ============================================================
     * 2. Load target, match k-mers, emit synteny blocks
     * ============================================================ */
    fprintf(stderr, "[msg] load target from %s\n", tarfile);

    FILE   *fp_ctg   = fopen(out_ctg, "w");
    gzFile  fp_block = gzopen(out_block, "w");
    if (!fp_ctg || !fp_block) {
        fprintf(stderr, "Cannot open output files\n"); return 1;
    }

    /* set of unique target kmers not in ref (matches Python's `kmer_new` set) */
    kh_strint_t *kmer_new = kh_strint_init();
    int    ctg_buf_size = 0, ctg_buf_cap = 1024;
    DARow *ctg_buf      = malloc(ctg_buf_cap * sizeof(DARow));
    char   prev_ctg[MAX_STR] = "";

    fp = gzopen(tarfile, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", tarfile); return 1; }

    while (gzgets(fp, line, sizeof(line))) {
        char ctg[MAX_STR], kmer[MAX_STR], dummy[MAX_STR], strand_c;
        long fro, to;
        if (sscanf(line, "%255s %ld %ld %255s %255s %c",
                   ctg, &fro, &to, kmer, dummy, &strand_c) != 6) continue;

        char *colon = strchr(kmer, ':');
        if (colon) *colon = '\0';

        khint_t k = kh_strmap_get(ref_dict, kmer);
        if (k == kh_end(ref_dict)) {
            int absent;
            khint_t kk = kh_strint_put(kmer_new, kmer, &absent);
            if (absent) kh_key(kmer_new, kk) = strdup(kmer);
            continue;
        }

        KmerEntry *e = (KmerEntry *)kh_val(ref_dict, k);
        e->count++;
        char strand = (strand_c == e->strand) ? '+' : '-';

        /* flush buffer when the contig changes */
        if (ctg_buf_size > 0 && strcmp(ctg, prev_ctg) != 0) {
            fprintf(stderr, "[msg] process contig %s\n", prev_ctg);
            SummaryRow *summ; int summ_n;
            BlockArr blocks;
            block_split(ctg_buf, ctg_buf_size, gap_skip_cut, gap_diff_cut,
                        &summ, &summ_n, &blocks);
            print_summary_line(fp_ctg, summ, summ_n);
            print_blocks(fp_block, &blocks);
            free(summ);
            free_blocks(&blocks);
            ctg_buf_size = 0;
        }
        copy_str(prev_ctg, ctg);

        if (ctg_buf_size == ctg_buf_cap) {
            ctg_buf_cap *= 2;
            ctg_buf      = realloc(ctg_buf, ctg_buf_cap * sizeof(DARow));
        }
        DARow *row = &ctg_buf[ctg_buf_size++];
        copy_str(row->ctg,     ctg);
        row->fro = fro; row->to = to;
        copy_str(row->ref_ctg, e->ctg);
        row->ref_fro = e->fro; row->ref_to = e->to;
        copy_str(row->kmer,    kmer);
        row->strand = strand;
    }
    gzclose(fp);

    /* flush the last contig */
    if (ctg_buf_size > 0) {
        fprintf(stderr, "[msg] process contig %s\n", ctg_buf[0].ctg);
        SummaryRow *summ; int summ_n;
        BlockArr blocks;
        block_split(ctg_buf, ctg_buf_size, gap_skip_cut, gap_diff_cut,
                    &summ, &summ_n, &blocks);
        print_summary_line(fp_ctg, summ, summ_n);
        print_blocks(fp_block, &blocks);
        free(summ);
        free_blocks(&blocks);
    }
    free(ctg_buf);
    fclose(fp_ctg);
    gzclose(fp_block);

    fprintf(stderr, "[msg] total %d new kmers found in target\n",
            (int)kh_size(kmer_new));
    for (khint_t ki = 0; ki < kh_end(kmer_new); ki++) {
        if (!kh_exist(kmer_new, ki)) continue;
        free((char *)kh_key(kmer_new, ki));
    }
    kh_strint_destroy(kmer_new);

    /* ============================================================
     * 3. Build CNV summary
     *    ref_dict.count now holds how many times each kmer was hit.
     * ============================================================ */
    kh_strmap_t *ctg_kmer_map = kh_strmap_init(); /* ctg -> CtgKmerList* */
    kh_intmap_t *cnv_map      = kh_intmap_init(); /* cp  -> kh_strint_t* */

    for (khint_t ki = 0; ki < kh_end(ref_dict); ki++) {
        if (!kh_exist(ref_dict, ki)) continue;
        KmerEntry *e = (KmerEntry *)kh_val(ref_dict, ki);

        /* --- ctg_kmer_map --- */
        khint_t k2 = kh_strmap_get(ctg_kmer_map, e->ctg);
        CtgKmerList *ckl;
        if (k2 == kh_end(ctg_kmer_map)) {
            ckl       = calloc(1, sizeof(CtgKmerList));
            copy_str(ckl->ctg, e->ctg);
            ckl->cap  = 16;
            ckl->data = malloc(ckl->cap * sizeof(KmerPos));
            int absent;
            k2 = kh_strmap_put(ctg_kmer_map, ckl->ctg, &absent);
            kh_val(ctg_kmer_map, k2) = ckl;
        } else {
            ckl = (CtgKmerList *)kh_val(ctg_kmer_map, k2);
        }
        if (ckl->size == ckl->cap) {
            ckl->cap *= 2;
            ckl->data = realloc(ckl->data, ckl->cap * sizeof(KmerPos));
        }
        KmerPos kp;
        kp.fro = e->fro; kp.to = e->to; kp.count = e->count;
        copy_str(kp.kmer, e->kmer);
        ckl->data[ckl->size++] = kp;

        /* --- cnv_map: cp -> (kh_strint_t: ctg -> n) --- */
        khint_t cp = (khint_t)e->count;
        khint_t k3 = kh_intmap_get(cnv_map, cp);
        kh_strint_t *inner;
        if (k3 == kh_end(cnv_map)) {
            inner = kh_strint_init();
            int absent;
            k3 = kh_intmap_put(cnv_map, cp, &absent);
            kh_val(cnv_map, k3) = inner;
        } else {
            inner = (kh_strint_t *)kh_val(cnv_map, k3);
        }
        int absent;
        khint_t k4 = kh_strint_put(inner, e->ctg, &absent);
        if (absent) kh_val(inner, k4) = 0;
        kh_val(inner, k4)++;
    }

    /* collect & sort cp values */
    int      n_cp  = (int)kh_size(cnv_map);
    khint_t *cp_arr = malloc(n_cp * sizeof(khint_t));
    {
        int i = 0;
        for (khint_t ki = 0; ki < kh_end(cnv_map); ki++)
            if (kh_exist(cnv_map, ki)) cp_arr[i++] = kh_key(cnv_map, ki);
    }
    qsort(cp_arr, n_cp, sizeof(khint_t), cmp_int);

    /* Drop any ctgs from ctg_ord whose kmers were all duplicate-removed
     * (mirrors Python: ref_ctg_set is built only from surviving kmers). */
    {
        int n_valid = 0;
        for (int i = 0; i < n_ctg_ord; i++) {
            if (kh_strmap_get(ctg_kmer_map, ctg_ord[i]) != kh_end(ctg_kmer_map))
                ctg_ord[n_valid++] = ctg_ord[i];
            else
                free(ctg_ord[i]);
        }
        n_ctg_ord = n_valid;
    }

    /* For kmer-summary.tsv, sort ref-ctgs alphabetically (matches Python's
     * sorted(ref_ctg_set)). The deletion-block pass below keeps ctg_ord in
     * first-appearance order. */
    char **ref_ctgs_sorted = malloc(n_ctg_ord * sizeof(char *));
    memcpy(ref_ctgs_sorted, ctg_ord, n_ctg_ord * sizeof(char *));
    qsort(ref_ctgs_sorted, n_ctg_ord, sizeof(char *), cmp_strptr);

    /* write kmer-summary.tsv */
    FILE *fp_kmer = fopen(out_kmer, "w");
    fprintf(fp_kmer, "Chrom");
    for (int i = 0; i < n_cp; i++) fprintf(fp_kmer, "\tcp%u", cp_arr[i]);
    fprintf(fp_kmer, "\n");

    for (int ri = 0; ri < n_ctg_ord; ri++) {
        fprintf(fp_kmer, "%s", ref_ctgs_sorted[ri]);
        for (int ci = 0; ci < n_cp; ci++) {
            khint_t k3 = kh_intmap_get(cnv_map, cp_arr[ci]);
            int val = 0;
            if (k3 != kh_end(cnv_map)) {
                kh_strint_t *inner = (kh_strint_t *)kh_val(cnv_map, k3);
                khint_t k4 = kh_strint_get(inner, ref_ctgs_sorted[ri]);
                if (k4 != kh_end(inner)) val = kh_val(inner, k4);
            }
            fprintf(fp_kmer, "\t%d", val);
        }
        fprintf(fp_kmer, "\n");
    }
    fclose(fp_kmer);
    free(ref_ctgs_sorted);
    free(cp_arr);

    /* ============================================================
     * 4. Deletion block detection
     *    For each ref-ctg: sort k-mers by fro, RLE by count.
     *    Output blocks where count == 0.
     * ============================================================ */
    int       cnv_n = 0, cnv_cap = 64;
    CnvBlock *cnv_blocks = malloc(cnv_cap * sizeof(CnvBlock));

    for (int ri = 0; ri < n_ctg_ord; ri++) {
        khint_t ki = kh_strmap_get(ctg_kmer_map, ctg_ord[ri]);
        if (ki == kh_end(ctg_kmer_map)) continue;
        CtgKmerList *ckl = (CtgKmerList *)kh_val(ctg_kmer_map, ki);

        qsort(ckl->data, ckl->size, sizeof(KmerPos), cmp_kmer_pos);

        int old_cnt = ckl->data[0].count;
        CnvBlock cur;
        copy_str(cur.ctg,        ckl->ctg);
        cur.fro   = ckl->data[0].fro;
        cur.to    = ckl->data[0].to;
        copy_str(cur.kmer_first, ckl->data[0].kmer);
        copy_str(cur.kmer_last,  ckl->data[0].kmer);
        cur.count = old_cnt;
        cur.n     = 1;

        for (int i = 1; i < ckl->size; i++) {
            int cnt = ckl->data[i].count;
            if (cnt == old_cnt) {
                cur.to = ckl->data[i].to;
                copy_str(cur.kmer_last, ckl->data[i].kmer);
                cur.n++;
            } else {
                if (cnv_n == cnv_cap) {
                    cnv_cap   *= 2;
                    cnv_blocks = realloc(cnv_blocks, cnv_cap * sizeof(CnvBlock));
                }
                cnv_blocks[cnv_n++] = cur;
                copy_str(cur.ctg,        ckl->ctg);
                cur.fro   = ckl->data[i].fro;
                cur.to    = ckl->data[i].to;
                copy_str(cur.kmer_first, ckl->data[i].kmer);
                copy_str(cur.kmer_last,  ckl->data[i].kmer);
                cur.count = cnt;
                cur.n     = 1;
                old_cnt   = cnt;
            }
        }
        if (cnv_n == cnv_cap) {
            cnv_cap   *= 2;
            cnv_blocks = realloc(cnv_blocks, cnv_cap * sizeof(CnvBlock));
        }
        cnv_blocks[cnv_n++] = cur;

        free(ckl->data); free(ckl);
    }
    kh_strmap_destroy(ctg_kmer_map);

    /*
     * Write all CNV blocks.
     * Columns: ctg  fro  to  kmers(first,last)  count  n
     */
    gzFile fp_del = gzopen(out_del, "w");
    if (!fp_del) { fprintf(stderr, "Cannot open %s\n", out_del); return 1; }
    for (int i = 0; i < cnv_n; i++) {
        CnvBlock *row = &cnv_blocks[i];
        char kmers[MAX_STR*2+2];
        snprintf(kmers, sizeof(kmers), "%s,%s", row->kmer_first, row->kmer_last);
        gzprintf(fp_del, "%s\t%ld\t%ld\t%s\t%d\t%d\n",
                 row->ctg, row->fro, row->to, kmers, row->count, row->n);
    }
    gzclose(fp_del);
    free(cnv_blocks);
    for (int i = 0; i < n_ctg_ord; i++) free(ctg_ord[i]);
    free(ctg_ord);

    /* ============================================================
     * Cleanup
     * ============================================================ */
    for (khint_t ki = 0; ki < kh_end(cnv_map); ki++) {
        if (!kh_exist(cnv_map, ki)) continue;
        kh_strint_destroy((kh_strint_t *)kh_val(cnv_map, ki));
    }
    kh_intmap_destroy(cnv_map);

    for (khint_t ki = 0; ki < kh_end(ref_dict); ki++) {
        if (!kh_exist(ref_dict, ki)) continue;
        free(kh_val(ref_dict, ki));
    }
    kh_strmap_destroy(ref_dict);

    fprintf(stderr, "[msg] done!\n");
    return 0;
}
