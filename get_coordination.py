from pylada.crystal import read, write, neighbors
import pylada.periodic_table as pt
from collections import defaultdict
import json
import pickle
import os
import numpy as np

d = json.load(open("data.json", 'r'))
cord = {}
coulomb = {"ZiZj/d": {}, 
           "1/d": {}}
for item in d:
    icsdstr = "{0:06d}".format(int(item["icsdnum"]))
    atoms = read.icsd_cif_a("/scratch/jyan/allcifs/icsd_%s.cif"%(icsdstr)) 
    ngh = defaultdict(float)
    uniqueatom = defaultdict(int)
    coul1 = defaultdict(float)
    coul2 = defaultdict(float)

    for i, atom in enumerate(atoms):
        ngh_n = neighbors(structure=atoms,nmax=1,center=atom.pos,tolerance=0.05)
        uniqueatom[atom.type] += 1
        ngh[atom.type] += len(ngh_n)
        
        Zi = getattr(pt, atom.type).atomic_number 
        # loop over neighbors
        for j in range(len(ngh_n)):
            Zj = getattr(pt, ngh_n[j][0].type).atomic_number 
            coul1[atom.type] += Zi * Zj / ngh_n[j][2]
            coul2[atom.type] += 1. / ngh_n[j][2]

    for atom in ngh.keys():
        ngh[atom] /= float(uniqueatom[atom])
        coul1[atom] /= float(uniqueatom[atom])
        coul2[atom] /= float(uniqueatom[atom])

    cord["%s"%(icsdstr)] = ngh
    coulomb["ZiZj/d"]["%s"%(icsdstr)] = coul1
    coulomb["1/d"]["%s"%(icsdstr)] = coul2
    print icsdstr

json.dump(cord, open('cord.json', 'w'))
json.dump(coulomb, open('coulomb.json', 'w'))
