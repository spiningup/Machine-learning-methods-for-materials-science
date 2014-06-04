import numpy as np
from filterdb import *
import shutil
import os
import commands

natoms = 100
nminatoms = 90

d = json.load(open('icsd.json', 'r'))
dnew = formula_with_multiple_spacegroups(d)    
d = filter_natoms(d, nmin=nminatoms, nmax=natoms)
for el in (['D', 'O', 'H', 'Re', 'Ta', 'Cs', 'Tl', 'Os', 'Tc']):
    d = filter_elements(d, el, incl=False)
print "number of atoms < %d and > %d"%(natoms, nminatoms), len(d)
#d = add_otherspacegroups(d, dnew)
#print "add other spacegropus", len(d)

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
    filename2 = "/scratch/jyan/ML_natoms_%s/icsd_%s.cif"%(natoms, key)

    if not os.path.isfile(filename2):
        shutil.copy2(filename, filename2)

#allcifs =  commands.getoutput("ls /scratch/jyan/ML_natoms_10/*.cif").split('\n')
#for cif in allcifs:
#    if cif[-10:-4] not in d.keys():
#        os.remove(cif) 
#        if os.path.isfile(cif[:-4]):
#            shutil.rmtree(cif[:-4])


