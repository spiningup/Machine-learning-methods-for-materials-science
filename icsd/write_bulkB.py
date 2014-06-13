import commands
import os
import json
import numpy as np
from pylada import vasp
from pylada.crystal import read

jobdir = os.getcwd()

a = os.listdir('.')
allcifs = [i for i in a if i.startswith("icsd_")]
print "number of calculated cifs", len(allcifs)

bulkBinput = []
for icsdno in allcifs:
    subdirs = os.listdir('%s/%s'%(jobdir, icsdno))
    E = []
    dirs = []
    for subdir in subdirs:
        a = vasp.Extract('%s/%s/%s'%(jobdir,icsdno, subdir))
        if a.success:
            E.append(a.total_energy.magnitude)
            dirs.append(subdir)
    if "non-magnetic" not in dirs: continue
    idx = np.argsort(E)[0]
    assert np.abs(E[idx] - min(E)) < 1e-8
    print "%s/%s/%s/CONTCAR"%(jobdir, icsdno, dirs[idx]), min(E)
    bulkBinput.append("%s/%s/%s/CONTCAR"%(jobdir, icsdno, dirs[idx]))
    
json.dump(bulkBinput, open("input_for_bulkB.json", "w"))

