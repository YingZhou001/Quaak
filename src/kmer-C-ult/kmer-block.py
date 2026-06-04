import sys,gzip
import numpy as np
from scipy import stats

if len(sys.argv) != 4:
    print("Usage: quick_block.py <ref.gz> <tar.gz> <outpref>", file=sys.stderr)
    sys.exit(1)


reffile = sys.argv[1]
tarfile = sys.argv[2]
outpref = sys.argv[3]


gap_skip_cut = 100000000
gap_diff_cut = 100000000


out_ctg_file = outpref + '.ctg-strand.tsv'
out_kmer_file = outpref + '.kmer-summary.tsv'
out_cnv_file = outpref + '.cnv.tsv.gz'
out_block_file = outpref + '.block.tsv.gz'

def cal_tau(x, y) :
    n = len(x)
    if n < 3 or n != len(y) : return np.nan, np.nan
    tau_statistic, p_value = stats.kendalltau(x, y)
    return tau_statistic, p_value


def block_split(DA, gap_skip_cut, gap_diff_cut) :

    #row in DA : [ctg, fro, to, ref_ctg, ref_fro, ref_to, kmer, strand]
    ret = []
    da_all = {}
    rgs = {}
    kmer_cnt = {}
    new_da_block = []
    for row in DA:
        ctg, fro, to, ref_ctg, ref_fro, ref_to, kmer, strand = row
        if kmer not in kmer_cnt : kmer_cnt[kmer] = 0
        kmer_cnt[kmer] += 1
        if ref_ctg not in da_all: da_all[ref_ctg] = [[],[],[]]
        da_all[ref_ctg][0].append(fro)
        da_all[ref_ctg][1].append(ref_fro)
        da_all[ref_ctg][2].append(strand)
        if ref_ctg not in rgs :
            rgs[ref_ctg] = [fro, to, ref_fro, ref_to]
        else :
            if fro < rgs[ref_ctg][0] : rgs[ref_ctg][0] = fro
            if to > rgs[ref_ctg][1] : rgs[ref_ctg][1] = to
            if ref_fro < rgs[ref_ctg][2] : rgs[ref_ctg][2] = ref_fro
            if ref_to > rgs[ref_ctg][3] : rgs[ref_ctg][3] = ref_to
    
    for row in DA:
        ctg, fro, to, ref_ctg, ref_fro, ref_to, kmer, strand = row
        if not new_da_block :
            new_da_block.append(row)
            continue
        row_old = new_da_block[-1]
        is_same_kmer_cnt = kmer_cnt[kmer] == kmer_cnt[row_old[6]]
        is_same_chr = row[3] == row_old[3]
        is_same_strand = row[7] == row_old[7]
        gap_tar = row[1] - row_old[1]
        gap_ref = row[4] - row_old[4]
        if gap_tar <= 0 :
            raise ValueError(f'Duplicated or unsorted rows in split_blocks: {row}')
        if gap_ref == 0 :
            is_same_order = False # force dup in different blocks
        else :
            is_same_order = gap_ref > 0 if row[7] == '+' else gap_ref < 0

        gap_ref_l = abs(gap_ref)
        # distance between kmers
        is_small_gap = gap_tar < gap_skip_cut and gap_ref_l < gap_skip_cut
        # distance difference between ref and tar
        is_same_gap = abs(gap_ref_l - gap_tar) < gap_diff_cut
        if is_same_chr and is_same_strand and is_same_order and is_small_gap and is_same_gap and is_same_kmer_cnt:
            new_da_block.append(row)
        else :
            ret.append(new_da_block)
            new_da_block = [row]
    if new_da_block: ret.append(new_da_block)


    summary = []
    for ref_ctg in da_all :
        x, y, strand = da_all[ref_ctg]
        coef, pvalue = cal_tau(x,y)
        n_plus = sum(1 for tmp in strand if tmp == '+')
        n_minus = len(x) - n_plus
        tar_fro, tar_to, ref_fro, ref_to = rgs[ref_ctg]
        summary.append([ctg, tar_fro, tar_to, ref_ctg, ref_fro, ref_to, coef, pvalue, n_plus, n_minus])

    return([summary, ret])


def print_summary_line(fp, summary) :
    for row in summary :
        outrow = [str(x) for x in row]
        print('\t'.join(outrow), file=fp)


def print_blocks(fp, blocks) :
    block_idx = 0
    for block in blocks:
        n = len(block)
        row0 = block[0]
        row1 = block[n-1]
        kmers = row0[6] + ','+ row1[6]
        tar_fro, tar_to = [min(row0[1], row1[1]), max(row0[2], row1[2])]
        ref_fro, ref_to = [min(row0[4], row1[4]), max(row0[5], row1[5])]
        ref_ctg = row0[3]
        n_plus = sum(1 for tmp in block if tmp[7] == '+')
        n_minus = n - n_plus
        ctg = row0[0]

        print(ctg, block_idx, tar_fro, tar_to, ref_ctg, ref_fro, ref_to, n_plus, n_minus, kmers, sep= '\t', file=fp)
        block_idx += 1

def print_cnv(fp, cnv, ref_ctg_set) :
    cp_arr = [cp for cp in cnv]
    cp_arr = sorted(cp_arr)

    head_str='Chrom\t' + '\t'.join(['cp'+str(x) for x in cp_arr])
    print(head_str, file=fp)
    for ref in sorted(ref_ctg_set) :
        smry = []
        for cp in cp_arr:
            if ref in cnv[cp] :
                smry.append(str(cnv[cp][ref]))
            else :
                smry.append('0')

        row_str = ref + '\t' + '\t'.join(smry)
        print(row_str, file=fp)


def print_cnv_blocks(fp, blocks) :
    for row in blcoks:
        print('\t'.join(row), file=fp)



print(f'[msg] load reference from {reffile}', file=sys.stderr)

ref_da_kmer_dict = {}
kmer_rm = set()

with gzip.open(reffile, 'rt') as fp:
    for line in fp:
        ctg, fro, to, kmer, _ , strand = line.split()
        if ':' in kmer : kmer = kmer.split(':')[0]
        if kmer not in ref_da_kmer_dict : 
            ref_da_kmer_dict[kmer] = [ctg, int(fro), int(to), strand, 0]
        else : kmer_rm.add(kmer)

for kmer in kmer_rm: del ref_da_kmer_dict[kmer]

print(f'[msg] {len(ref_da_kmer_dict)} kmer loaded, {len(kmer_rm)} removed due to duplication in ref', file=sys.stderr)
del kmer_rm



print(f'[msg] load target from {tarfile}', file=sys.stderr)


with open(out_ctg_file, 'wt') as fp_ctg, \
        gzip.open(out_block_file, 'wt') as fp_block:

    ctg_buf = []
    kmer_new = set()
    with gzip.open(tarfile, 'rt') as fp:
        for line in fp:
            ctg, fro, to, kmer, _, strand = line.split()
            if ':' in kmer : kmer = kmer.split(':')[0]
            if kmer not in ref_da_kmer_dict : 
                kmer_new.add(kmer)
                continue
            ref_ctg, ref_fro, ref_to, ref_strand, _ = ref_da_kmer_dict[kmer]
            ref_da_kmer_dict[kmer][4] += 1 
            strand = '+' if strand == ref_strand else '-'
            row = [ctg, int(fro), int(to), ref_ctg, ref_fro, ref_to, kmer, strand]
            if not ctg_buf :
                ctg_buf = [row]
                continue


            old_ctg = ctg_buf[-1][0]
            if old_ctg != ctg :
                # process ctg_buf
                print(f'[msg] process contig {old_ctg}', file=sys.stderr)
                summary, blocks = block_split(ctg_buf, gap_skip_cut, gap_diff_cut)
                print_summary_line(fp_ctg, summary)
                print_blocks(fp_block, blocks)
                ctg_buf = [row]
            else :
                ctg_buf.append(row)

    if ctg_buf :
        # process ctg_buf
        ctg = ctg_buf[0][0]
        print(f'[msg] process contig {ctg}', file=sys.stderr)
        summary, blocks = block_split(ctg_buf, gap_skip_cut, gap_diff_cut)
        print_summary_line(fp_ctg, summary)
        print_blocks(fp_block, blocks)

print(f'[msg] total {len(kmer_new)} new kmers not used in target', file=sys.stderr)

# process kmer CNV summary
cnv = {}
ref_da_ctg_dict = {}
ref_ctg_set = set()

for kmer in ref_da_kmer_dict :
    ctg, fro, to, _, cnt = ref_da_kmer_dict[kmer]
    if ctg not in ref_da_ctg_dict : ref_da_ctg_dict[ctg] = []
    ref_da_ctg_dict[ctg].append([fro, to, kmer, cnt])
    ref_ctg_set.add(ctg)
    if cnt not in cnv : cnv[cnt] = {}
    if ctg not in cnv[cnt]: cnv[cnt][ctg] = 0
    cnv[cnt][ctg] += 1

with open(out_kmer_file, 'wt') as fp_kmer :
    print_cnv(fp_kmer, cnv, ref_ctg_set)

# detecting deletion blocks


cnv_blocks = []
for ctg in ref_da_ctg_dict :
    tmp_da = sorted(ref_da_ctg_dict[ctg], key = lambda x: x[0])
    fro, to, kmer, old_cnt = tmp_da[0]
    tmp_arr = [ctg, fro, to, kmer, kmer, old_cnt, 0]
    for fro, to, kmer, cnt in tmp_da :
        if cnt == old_cnt :
            tmp_arr[2] = to
            tmp_arr[4] = kmer
            tmp_arr[6] += 1
        else :
            cnv_blocks.append(tmp_arr)
            tmp_arr = [ctg, fro, to, kmer, kmer, cnt, 1]
            old_cnt = cnt

    if tmp_arr: cnv_blocks.append(tmp_arr)


with open(out_cnv_file, 'wt') as fp_cnv :
    print_cnv_blocks(fp_cnv, cnv_blocks)

print(f'[msg] done!', file=sys.stderr)
