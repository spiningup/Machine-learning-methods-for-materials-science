from pylada.crystal import read
import json
import numpy as np

exptvol = {}

d = json.load(open("data.json", 'r'))
for item in d:
    natoms = int(item["numatom"])
    masses = np.array(item["atommasses_amu"])            
    names = item["atomnames"]
    name = names[np.argsort(masses)[0]]
    icsdstr = "{0:06d}".format(int(item["icsdnum"]))
    try:
        atoms = read.icsd_cif_a("/home/jyan/icsd/%s/icsd_%s.cif"%(name, icsdstr))
    except:
        continue
    vol = atoms.volume.magnitude / len(atoms)
    exptvol["%s"%(icsdstr)] = vol

json.dump(exptvol, open('exptvol.json', 'w'))
