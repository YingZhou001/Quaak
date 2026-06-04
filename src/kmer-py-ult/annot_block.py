import sys,gzip
import numpy as np
from collections import Counter

import mwis_permutation as mp

# load confident contigs
ctg_strand = {}

with open(sys.argv[1], 'rt') as fp:
    for line in fp:
        ctg, ref, tau = line.split()
        tau = float(tau)
        if ctg not in ctg_strand : ctg_strand[ctg] = {}
        ctg_strand[ctg][ref] = '+' if tau > 0 else '-' 


da = {}
with gzip.open(sys.argv[2], 'rt') as fp:
    for line in fp:
        if line[0] == '#' : continue
        llst = line.split()
        ctg = llst[0]
        ref = llst[4]
        if ctg in ctg_strand and ref in ctg_strand[ctg]:
            if ctg not in da: da[ctg] = {}
            if ref not in da[ctg]: da[ctg][ref] = []
            da[ctg][ref].append(llst)

ref_path = {}
with gzip.open(sys.argv[3], 'rt') as fp:
    for line in fp:
        if line[0] == '#' : continue
        llst = line.split()
        ref = llst[0]
        fro = int(llst[1])
        if ref not in ref_path :
            ref_path[ref] = []
        ref_path[ref].append(fro)

for ref in ref_path:
    ref_path[ref] = np.array(ref_path[ref])


# reverse the backward strand
for ctg in da:
    for ref in da[ctg]:
        strand = ctg_strand[ctg][ref]
        fro = da[ctg][ref][0][2]
        to = da[ctg][ref][-1][3]
        ctg0 = ctg + ':' + fro + '-' + to
        if strand == '-' :
            da[ctg][ref] = da[ctg][ref][::-1]
            for i in range(len(da[ctg][ref])) :
                row = da[ctg][ref][i]
                row[0] = ctg0
                ## reverse the position
                a,b = int(row[2]),int(row[3])
                row[2] = str(int(to) + 1 - b)
                row[3] = str(int(to) + 1 - a)
                ## reverse the strand direction
                row[7], row[8] = row[8], row[7]
                da[ctg][ref][i] = row


def find_anchor_blocks(blocks) :
    # use the gapless input, output the idx set
    # three conditions to define the anchor blocks
    ## 1) positive, longer than 20 kmers
    ## 2) not duplicated 
    ## 3) use boundary position to define the duplication
    n = len(blocks)
    ret_idx_set = set()
    left_arr = []
    right_arr = []
    for i in range(n) :
        n_pos = int(blocks[i][7])
        if n_pos > 20 : ret_idx_set.add(i)
        left_arr.append(blocks[i][5])
        right_arr.append(blocks[i][6])

    left_counts = Counter(left_arr)
    right_counts = Counter(right_arr)
    dups_set = set()
    for item, count in left_counts.items():
        if count > 1 : dups_set.add(item)

    for item, count in right_counts.items():
        if count > 1 : dups_set.add(item)

    for i in range(n) :
        c1 = left_arr[i] in dups_set
        c2 = right_arr[i] in dups_set
        if c1 or c2 : ret_idx_set.discard(i)

    ret_idx_arr = sorted(ret_idx_set)

    # remove crossing blocks
    # think about the problem here
    # we have two arr x and y
    # x = [1, 2, 3, 4, 5, 6, 7], and y = [1, 2, 4, 3, 5, 6, 7] 
    # which is a permutation of x
    # with weight block length for each item
    # we need to find out the minimal cost to items so that the
    # remain two arrays are both in strictly increasing order

    m = len(ret_idx_arr) 
    x = []
    xy = []
    w = []
    for k in range(m) :
        i = ret_idx_arr[k]
        x.append(int(blocks[i][3]))
        xy.append([int(blocks[i][3]), int(blocks[i][5])])
        w.append(int(blocks[i][3]) - int(blocks[i][2]))
    y = [z[0] for z in sorted(xy, key = lambda x: x[1])]
    if not x : return []
    best, kept, removed = mp.max_weight_kept(x, y, w)
    ret_idx_arr_final = [ret_idx_arr[i] for i in kept]

    return(ret_idx_arr_final)


# find out anchor blocks
anchor_blocks = {}
for ctg in da:
    if ctg not in anchor_blocks : anchor_blocks[ctg] = {}
    for ref in da[ctg]:
        anchor_blocks[ctg][ref] = find_anchor_blocks(da[ctg][ref])



def search_anchor_boundaries(i, dat, A_idx_set) :
    n = len(dat)
    row = dat[i]
    left = -1
    right = n
    if i in A_idx_set :
        left = i
        right = i
    else :
        if i > 0 :
            left = i - 1
            while left >= 0 :
                if left in A_idx_set : break
                else : left -= 1
        
        if i < n-1 :
            right = i + 1
            while right < n :
                if right in A_idx_set: break
                else : right += 1

    if left == -1 : ref_left_end = 'NA'
    else : ref_left_end = dat[left][6]
    if right == n : ref_right_end = 'NA'
    else : ref_right_end = dat[right][5]
    if left == right :
        ref_left_end = dat[right][5]
        ref_right_end = dat[left][6]

    return([left, right, ref_left_end, ref_right_end])


def find_dup_set(dat) :
    ref_pos_arr = []
    for row in dat :
        ref_pos_arr.append(row[5])
        ref_pos_arr.append(row[6])
    
    counts = Counter(ref_pos_arr)
    dups = [item for item, count in counts.items() if count > 1]
    return(dups)


def restore_inverted_ctg(row) :
    ctg, rg0 = row[0].split(':')
    f0,t0 = [int(x) for x in rg0.split('-')]
    f1,t1 = [int(row[2]), int(row[3])]
    row[0] = ctg
    row[2] = str(t0 + 1 - t1)
    row[3] = str(t0 + 1 - f1)
    return(row)


def find_del_boundary(buf, path) :
    n = len(path)
    in_a = buf[0]
    in_b = buf[-1]
    in_left = str(path[in_a])
    in_right = str(path[in_b] + 31)
    out_a = in_a - 1
    out_b = in_b + 1
    # DEL drop at both ends
    if out_a < 0 : return []
    else : out_left = str(path[out_a])
    if out_b > n-1 : return []
    else : out_right = str(path[out_b] + 31)
    in_kmer_num = str(len(buf))
    return([in_left, in_right, out_left, out_right, in_kmer_num])




for ctg in da:
    for ref in da[ctg]:
        out0 = [] # well-define region
        X = da[ctg][ref]
        A_idx_set = set(anchor_blocks[ctg][ref])
        ref_cnv_path = np.zeros(len(ref_path[ref]))
        dups = find_dup_set(X)
        n = len(X)
        for i in range(n):
            row = X[i]
            tags = []
            l_idx, r_idx, rl_end, rr_end = search_anchor_boundaries(i, X, A_idx_set)
            # check anchor block
            if i in A_idx_set : tags.append('ANC')
            # check insertion
            cl = rl_end != 'NA' and int(row[6]) < int(rl_end)
            cr = rr_end != 'NA' and int(row[5]) > int(rr_end)
            # check insertion
            if cl or cr : tags.append('INS')
            # check inversion
            if int(row[8]) > 0 : tags.append('INV')
            # check duplication
            if row[5] in dups or row[6] in dups : tags.append('DUP')

            if not tags: tags = ['other']

            if ':' in row[0] :
                row = restore_inverted_ctg(row)

            ctg_in_fro = row[2]
            ctg_in_to = row[3]
            ref_in_fro = row[5]
            ref_in_to = row[6]
            ref_out_fro = rl_end
            ref_out_to = rr_end
            in_kmer_num = str(int(row[7]) + int(row[8]))
            outstr = [row[0], ctg_in_fro, ctg_in_to, 
                      ref, ref_in_fro, ref_in_to, ref_out_fro, ref_out_to,
                      in_kmer_num, ','.join(tags)]
            r_fro = int(row[5])
            r_to = int(row[6])
            o = (r_fro <= ref_path[ref]) & (ref_path[ref] < r_to)
            ref_cnv_path[o] += 1

            out0.append(outstr)
        
        buf = []
        rn = len(ref_cnv_path)
        for i in range(rn) :
            if ref_cnv_path[i] == 0 :
                buf.append(i)
            else :
                if buf :
                    ret = find_del_boundary(buf, ref_path[ref])
                    if ret : 
                        ref_in_fro, ref_in_to, ref_out_fro, ref_out_to, in_kmer_num = ret
                        outstr = [ctg, 'NA', 'NA', ref, 
                                  ref_in_fro, ref_in_to,
                                  ref_out_fro, ref_out_to,
                                  in_kmer_num, 'DEL']
                        out0.append(outstr)
                    buf = []
        
        if buf :
            ret = find_del_boundary(buf, ref_path[ref])
            if ret : 
                ref_in_fro, ref_in_to, ref_out_fro, ref_out_to, in_kmer_num = ret
                outstr = [ctg, 'NA', 'NA', ref, 
                          ref_in_fro, ref_in_to,
                          ref_out_fro, ref_out_to,
                          in_kmer_num, 'DEL']
                out0.append(outstr)


        if out0: 
            #out0 = sorted(out0, key = lambda x: int(x[3].split('..')[0]))
            for row in out0:
                print('\t'.join(row))
