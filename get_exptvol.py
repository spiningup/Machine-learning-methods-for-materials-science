from pylada.crystal import read, write
import json
import os
import numpy as np
import shutil
from collections import Counter


exptvol = {}

d = json.load(open("data.json", 'r'))
for item in d:
    natoms = int(item["numatom"])
    masses = np.array(item["atommasses_amu"])            
    names = item["atomnames"]
    name = names[np.argsort(masses)[0]]
    icsdstr = "{0:06d}".format(int(item["icsdnum"]))
    distinctnames = Counter(names)

    if not os.path.isfile("/home/jyan/icsd/%s/icsd_%s.cif"%(name, icsdstr)):
        for name in distinctnames.keys():
            if os.path.isfile("/home/jyan/icsd/%s/icsd_%s.cif"%(name, icsdstr)):
                break
            

#    print name, icsdstr
    atoms = read.icsd_cif_a("/home/jyan/icsd/%s/icsd_%s.cif"%(name, icsdstr))        
    shutil.copy2("/home/jyan/icsd/%s/icsd_%s.cif"%(name, icsdstr), 'structures/icsd_%s.cif'%(icsdstr))

    vol = atoms.volume.magnitude / len(atoms)
    exptvol["%s"%(icsdstr)] = vol

json.dump(exptvol, open('exptvol.json', 'w'))
