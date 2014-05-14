import numpy as np
import json
import os
import commands
from pylada.crystal import read

d = json.load(open('icsd.json', 'r'))

dnew = {}
ii = 0
for key, item in d.items():

    formula = item[0].replace("'",'')
    for i in formula:
        if i.isdigit():
            formula = formula.replace(i, '')
    formula = formula.split()

    for el in formula:
        if os.path.isfile("/home/jyan/icsd/%s/icsd_%s.cif"%(el, key)):
            break
    filename = "/home/jyan/icsd/%s/icsd_%s.cif"%(el, key)

    try:
        atoms = read.icsd_cif_a(filename)
    except:
        print key, item
        continue

    dnew[key] = item + [len(atoms), ]
    print ii, dnew[key]
    ii += 1

json.dump(dnew, open("icsd_wi_natoms.json", 'w'))

