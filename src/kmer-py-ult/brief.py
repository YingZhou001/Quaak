import sys,gzip


oneline = []

print("CC\tComment lines")

print("CC\tKS\tcategory,autosome,auto_percentage, chrX,chrY")
print("CC\tIC\tqctg,qfro,qto,rctg,rfro,rto,tau,pvalue,n_kmer_pos,n_kmer_neg")
print("CC\tBR\trctg,fro,to,boundary_kmer,cnv,nkmer")
print("CC\tLS\tqctg,qfro,qto,rctg,rin_fro,rin_to,rout_fro,rout_to,n_kmer,type")
print("CC\tPS\tcategory,count")
print("CC\tSN\ttag,count,annotation")

print("CC")
print("CC")

print("CC\t#Summary of unique kmer count")


out_A = []
out_X = []
out_Y = []
with open(sys.argv[1], 'rt') as fp:
    for line in fp:
        llst = line.split()
        if llst[0] == 'Chrom' :
            pass
        else :
            if 'X' in llst[0] :
                out_X = [int(x) for x in llst[1:]]
            elif 'Y' in llst[0] :
                out_Y = [int(x) for x in llst[1:]]
            else :
                if not out_A :
                    out_A = [int(x) for x in llst[1:]]
                else :
                    for i in range(len(out_A)) :
                        out_A[i] += int(llst[i+1])


A_missing = out_A[0]
A_match = out_A[1]
A_dup = sum(out_A[2:])
A_tot = sum(out_A)


X_missing = out_X[0]
X_match = out_X[1]
X_dup = sum(out_X[2:])
X_tot = sum(out_X)


Y_missing = out_Y[0]
Y_match = out_Y[1]
Y_dup = sum(out_Y[2:])
Y_tot = sum(out_Y)

missing = A_missing + X_missing + Y_missing
match = A_match + X_match + Y_match
dup = A_dup + X_dup + Y_dup
tot = A_tot + X_tot + Y_tot


print(f'KS\t[M]issing\t{A_missing}\t{100.0*A_missing/A_tot:.02f}%\t{X_missing}\t{Y_missing}')
print(f'KS\t[S]ingle-copy\t{A_match}\t{100.0*A_match/A_tot:.02f}%\t{X_match}\t{Y_match}')
print(f'KS\t[D]uplicated\t{A_dup}\t{100.0*A_dup/A_tot:.02f}%\t{X_dup}\t{Y_dup}')
print(f'KS\t[T]otal-count\t{A_tot}\t100%\t{X_tot}\t{Y_tot}')
print("CC")

oneline.append(["KMA", A_missing, "missing kmers in autosome"])
oneline.append(["KMX", X_missing, "missing kmers in chrX"])
oneline.append(["KMY", Y_missing, "missing kmers in chrY"])
oneline.append(["KSA", A_match, "single-copy kmers in autosome"])
oneline.append(["KSX", X_match, "single-copy kmers in chrX"])
oneline.append(["KSY", Y_match, "single-copy kmers in chrY"])
oneline.append(["KDA", A_dup, "duplicated kmers in autosome"])
oneline.append(["KDX", X_dup, "duplicated kmers in chrX"])
oneline.append(["KDY", Y_dup, "duplicated kmers in chrY"])

print("CC")

def count_file_lines(filepath) :
    with open(filepath, 'r') as f:
        count = 0
        for count, line in enumerate(f, 1): pass

    return(count)

taucut = sys.argv[2]
conf_ctgs_count = count_file_lines(sys.argv[3])
susp_ctgs_count = count_file_lines(sys.argv[4])
resc_ctgs_count = count_file_lines(sys.argv[5])

print(f"CC\t#Determination of contig-reference pair strand, taucut = {taucut}")
print(f"PS\tConfident-pair\t{conf_ctgs_count}")
print(f"PS\tLow-confi-pair\t{susp_ctgs_count}")
print(f"PS\tManual-pair\t{resc_ctgs_count}")
print("CC")

oneline.append(["CCS", conf_ctgs_count, "Contig-ref pair with confident strand direction"])
oneline.append(["CLS", susp_ctgs_count, "contig-ref pair with low-confident strand direction"])
oneline.append(["CMS", resc_ctgs_count, "contig-ref pair with maunually adjusted strand direction"])


print("CC")


l_min, k_min = [int(x) for x in sys.argv[6].split(',')]


print(f"CC\t#Summary of large SV events: len_min = {l_min}bp, kmer_min = {k_min}")


print("CC\t##Inter-chromosome events, filtered by the number of kmers within inner boundary")

ctg_anc = {}
inter_chr_ctg = set()
with open(sys.argv[7], 'rt') as fp:
    for line in fp:
        llst = line.split()
        ctg = llst[0]
        ref = llst[3]
        l = int(llst[2]) - int(llst[1])
        nkmer = int(llst[8]) + int(llst[9])
        if ctg not in ctg_anc : ctg_anc[ctg] = []
        if nkmer > k_min :
            ctg_anc[ctg].append(llst)

out = []
for ctg in ctg_anc :
    if len(ctg_anc[ctg]) > 1 :
        inter_chr_ctg.add(ctg)
        for row in ctg_anc[ctg] :
            out.append('\t'.join(row))


print(f"CC\t{len(inter_chr_ctg)} contigs found with large Inter-chromosome events")

if out :
    for row in out: print('IC\t'+row)

print("CC")

oneline.append(["CIE", len(inter_chr_ctg), "contigs with large Inter-chromosome events"])


print("CC")



print("CC\t##CNV blocks, filtered by inner boundary")

out = []
n_mis_blk = 0
n_dup_blk = 0
with gzip.open(sys.argv[8], 'rt') as fp:
    for line in fp:
        llst = line.split()
        ctg = llst[0]
        fro = llst[1]
        to = llst[2]
        cp = llst[4]
        nkmer = llst[5]
        if cp!= '1' and int(to) - int(fro) >= l_min and int(nkmer) >= k_min :
            out.append(llst)
            if cp == '0' : n_mis_blk += 1
            else : n_dup_blk += 1

print(f"CC\t{n_mis_blk} missing and {n_dup_blk} duplication blocks in the target genome")
if out:
    for row in out: print("BR\t" + '\t'.join(row))

oneline.append(["BMR", n_mis_blk, "Missing blocks in reference genome"])
oneline.append(["BDR", n_dup_blk, "Duplicated blocks in reference genome"])
print("CC")


da = {}
sel_region = {}
with gzip.open(sys.argv[9], 'rt') as fp:
    for line in fp:
        llst = line.split()
        ctg = llst[0]
        ref = llst[3]
        fro = llst[6]
        to = llst[7]
        if fro == 'NA' and to == 'NA' : continue
        rg = fro + '..' + to
        if fro == 'NA' or to == 'NA' :
            sv_l = l_min + 1 
        else :
            sv_l = int(to) - int(fro)
        tag = llst[9]
        if tag != 'ANC' and tag != 'DEL' and sv_l > l_min :
            if ctg not in sel_region : sel_region[ctg] = {}
            if ref not in sel_region[ctg]: sel_region[ctg][ref] = {}
            if rg not in sel_region[ctg][ref] : 
                sel_region[ctg][ref][rg] = [[], 0, []]
            sel_region[ctg][ref][rg][1] += int(llst[8])
            sel_region[ctg][ref][rg][2].append(llst[9])
            if 'NA' in rg :
                sel_region[ctg][ref][rg][0].append(llst[1])
                sel_region[ctg][ref][rg][0].append(llst[2])
            else :
                sel_region[ctg][ref][rg][0] = [sv_l]

        if ref not in da: da[ref] = []
        da[ref].append(llst)


n_large_sv = 0
n_simple_sv = 0
n_complex_sv = 0

for ctg in sel_region:
    for ref in sel_region[ctg]:
        for rg in sel_region[ctg][ref]:
            if 'NA' in rg :
                tmp = []
                for x in sel_region[ctg][ref][rg][0]:
                    if x != 'NA' :
                        tmp.append(int(x))
                sv_l = max(tmp) - min(tmp)
                if sv_l < l_min : 
                    sel_region[ctg][ref][rg] = []
                    continue
            n_kmer = sel_region[ctg][ref][rg][1]
            if n_kmer < k_min :
                sel_region[ctg][ref][rg] = []
                continue
            n_large_sv += 1
            if len(sel_region[ctg][ref][rg][2]) > 1 : n_complex_sv += 1
            else : n_simple_sv += 1

#print(sel_region)

print("CC")
print("CC\t##Large SVs: (from strand confident ctgs + recused if any)")
print("CC\tfiltered by outer boundary due to complexity, check *.sel-sv.pdf for visualization")
print("CC\tif any one of outer boundary not defined, inner range on the target contig is used. check *.sel-sv.pdf for visualization")
print("CC\tDUP are within contig duplication")
print(f"CC\t{n_large_sv} large SVs are found: {n_simple_sv} simple SVs and {n_complex_sv} complex SVs")


oneline.append(["SVS", n_simple_sv, "large simple SVs"])
oneline.append(["SVC", n_complex_sv, "large complex SVs"])


out = []

for ref in da:
    for row in da[ref] :
        ctg = row[0]
        rg = row[6] + '..' + row[7]
        if ctg in sel_region and ref in sel_region[ctg] \
           and rg in sel_region[ctg][ref] \
           and sel_region[ctg][ref][rg] :
            out.append(row)

if out:
    for row in out:
        print("LS\t" + '\t'.join(row))



print("CC")
print("CC\t#Summary numbers")
for row in oneline :
    print("SN\t" + '\t'.join([str(x) for x in row]))
