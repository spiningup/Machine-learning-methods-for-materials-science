import numpy as np
from filterdb import *
import shutil

d = json.load(open('icsd.json', 'r'))
dnew = formula_with_multiple_spacegroups(d)    
d = filter_natoms(d, nmax=10)
for el in (['D', 'O', 'H']):
    d = filter_elements(d, el, incl=False)
print "number of atoms < 10", len(d)
d = add_otherspacegroups(d, dnew)
print "add other spacegropus", len(d)

dnonmag, dmag = filter_magnetic_elements(d)
print "total, magnetic, non-magnetic", len(d), len(dmag), len(dnonmag)

# print number of magnetic atoms for each formula
#for key, item in dmag.items(): 
#    print key, item[0], item[3]

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
    shutil.copy2(filename, "/scratch/jyan/ML_natoms_10/icsd_%s.cif"%(key))
    
