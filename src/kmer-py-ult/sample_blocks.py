import sys

# load contigs

ctgs = set()

with open(sys.argv[1], 'rt') as fp:
    for line in fp:
        ctg,ref,tau = line.split()
        ctgs.add(ctg+','+ref)

for line in sys.stdin:
    llst = line.split()
    k = llst[0] + ',' + llst[4]
    if k in ctgs :
        print(line, end = '')
