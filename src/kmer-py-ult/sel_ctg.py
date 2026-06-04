import sys,gzip


tau_cut = float(sys.argv[1])
cfile = sys.argv[2]
sfile = sys.argv[3]

out = {}

fp_c = open(cfile, 'wt')
fp_s = open(sfile, 'wt')

for line in sys.stdin:
    llst = line.split()
    ctg = llst[0]
    ref = llst[3]
    tau = llst[6]
    if llst[8] == '0' or llst[9] == '0' :
        if llst[8] == '0' :
            print(ctg, ref, -1, sep=' ', file = fp_c)
        else :
            print(ctg, ref, 1, sep=' ', file = fp_c)
    elif tau != 'nan' :
        tau = float(tau)
        if abs(tau) > tau_cut :
            print(ctg, ref, tau, sep=' ', file = fp_c)
        else :
            print(ctg, ref, tau, sep=' ', file = fp_s)

fp_c.close()
fp_s.close()
