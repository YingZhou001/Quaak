import sys


gid = sys.argv[1]
ctg_strand_file = sys.argv[2]


REF,rgs = gid.split(':')
FRO,TO = [x for x in rgs.split('-')]
#if FRO == 'NA' : FRO = 0
#else : FRO = int(FRO)
#if TO == 'NA' : TO = 9e9
#else : TO = int(TO)
#l = TO - FRO

ctg_strand = {}
with open(ctg_strand_file, 'rt') as fp:
    for line in fp:
        ctg,ref,tau = line.split()
        if ref != REF : continue
        if float(tau) > 0 : ctg_strand[ctg] = '+'
        else : ctg_strand[ctg] = '-'



sel_ctgs = set()
da = {}

for line in sys.stdin:
    llst = line.split()
    ctg = llst[0]
    r = llst[3]
    if ctg not in ctg_strand: continue
    if r != REF : continue
    if ctg not in da: da[ctg] = []
    da[ctg].append(llst)

# cut the block

def is_overlap(a, b) :
    if a[0] > b[1] or a[1] < b[0] : return False
    return True

out = []
for ctg in da:
    da_x = da[ctg]
    n = len(da_x)
    sel_idx = set()
    for i in range(n):
        if da_x[i][6] == FRO and da_x[i][7] == TO \
           or da_x[i][6] == TO or da_x[i][7] == FRO :
            sel_idx.add(i)
    if not sel_idx: continue
    low, up = min(sel_idx), max(sel_idx) + 1
    for row in da_x[low:up] :
        rowstr = '\t'.join(row)
        if rowstr not in out: out.append(rowstr)

for rowstr in out:
    print(rowstr)
