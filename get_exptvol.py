from pylada.crystal import read, write
import json
import pickle
import os
import numpy as np
import shutil
from collections import Counter
import commands

exptvol = {}

d = json.load(open("data.json", 'r'))
for item in d:
    natoms = int(item["numatom"])
    masses = np.array(item["atommasses_amu"])            
    names = item["atomnames"]
    name = names[np.argsort(masses)[0]]
    icsdstr = "{0:06d}".format(int(item["icsdnum"]))
    distinctnames = Counter(names)

    # get volume per atom
    atoms = read.icsd_cif_a("/scratch/jyan/allcifs/icsd_%s.cif"%(icsdstr)) 
    vol = atoms.volume.magnitude / len(atoms)

    # get a,b,c, alpha,beta,gamma
    values = np.zeros(7)
    for i, name in enumerate(["_cell_length_a", "_cell_length_b", "_cell_length_c", 
                             "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma", "_cell_volume"]):
        values[i] = float(commands.getoutput("grep %s /scratch/jyan/allcifs/icsd_%s.cif"%(name, icsdstr)).split()[1].split("(")[0])

    values[0:3] /= (values[6] / vol)**(1./3.)
    values[6] = vol

    exptvol["%s"%(icsdstr)] = values
    print icsdstr, values

pickle.dump(exptvol, open('exptvol.pkl', 'wb'), -1)
#json.dump(exptvol, open('exptvol.json', 'w'))
