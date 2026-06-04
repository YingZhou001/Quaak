import sys,gzip


l_min, k_min = [int(x) for x in sys.argv[1].split(',')]

da = {}
sel_region = {}
with gzip.open(sys.argv[2], 'rt') as fp:
    for line in fp:
        llst = line.split()
        ctg = llst[0]
        ref = llst[3]
        rg = llst[6] + '-' + llst[7]
        #if 'NA' in rg: continue
        #fro,to = [int(x) for x in llst[6:8]]
        tag = llst[9]
        #if tag != 'ANC' and to - fro > l_min :
        if tag not in ['ANC', 'DEL'] :
            if ctg not in sel_region: sel_region[ctg] = {}
            if ref not in sel_region[ctg]: sel_region[ctg][ref] = {}
            if rg not in sel_region[ctg][ref] : 
                sel_region[ctg][ref][rg] = [-1, -1, 0]
            if llst[6] != 'NA' : fro = int(llst[6])
            elif llst[4] != 'NA' : fro = int(llst[4])
            else : continue
            if llst[7] != 'NA' : to = int(llst[7])
            elif llst[5] != 'NA' : to = int(llst[5])
            else : continue
            a = sel_region[ctg][ref][rg][0]
            b = sel_region[ctg][ref][rg][1]
            if a < 0 :
                sel_region[ctg][ref][rg][0] = fro
            else :
                sel_region[ctg][ref][rg][0] = min(fro, a)
            if b < 0 :
                sel_region[ctg][ref][rg][1] = to
            else :
                sel_region[ctg][ref][rg][1] = max(to, b)

            sel_region[ctg][ref][rg][2] += int(llst[8])

        if ref not in da: da[ref] = []
        da[ref].append(llst)

out = []

for ref in da:
    for row in da[ref] :
        ctg = row[0]
        rg = row[6] + '-' + row[7]
        if ctg in sel_region \
           and ref in sel_region[ctg] \
           and rg in sel_region[ctg][ref] :
            a,b,nk = sel_region[ctg][ref][rg]
            if b - a > l_min and nk > k_min :
                print(ref + ':' + rg, ctg, sep='\t')
