import numpy as np
from collections import defaultdict
from pylab import *
import json

f = open('bandgap_expt.dat', 'r')
items = f.readlines()
f.close()

bandgaps = defaultdict(list)
for item in items:
    index, A, i, B, j, gap, T = item.split()
    formula = " ".join([A, i, B, j])
    formula = formula.replace(" NA","")
    if gap != "NA":
        bandgaps[formula].append([int(index), float(gap), float(T)])

exptgap = []
for key, value in bandgaps.items():
    a = np.array(value)
    minT = np.argsort(a[:,2])[0] # minimum T
    gap = a[minT, 1]
    exptgap.append({"formula": key, "gap": gap})

json.dump(exptgap, open('exptgap.json', 'w'))
    
