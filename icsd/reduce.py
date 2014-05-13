import numpy as np
import os
import commands
import json

dir0 = '/home/jyan/icsd'
dirs = os.listdir(dir0)

excl_els = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
            "Tm", "Yb", "Lu", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", 
            "Cf", "Es", "Fm", "Md", "No", "Lr",
            "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
            "He", "Ne", "Ar", "Kr", "Xe", "Rn", "Po", "At", "Fr", "Ra"]

icsdat = {}
existset = []

for dir in dirs:
    print dir
    icsdfiles =  os.listdir('%s/%s'%(dir0,dir))
    for icsd in icsdfiles:
        try:
            formula = commands.getoutput('grep _chemical_formula_sum %s/%s/%s'%(dir0, dir, icsd)).split()[1:]
            formula = " ".join(sorted(formula))
            vol =  commands.getoutput('grep _cell_volume %s/%s/%s'%(dir0, dir, icsd)).split()[1]
            groupnum =  commands.getoutput('grep _symmetry_Int_Tables_number %s/%s/%s'%(dir0, dir, icsd)).split()[1]
        except:
            continue

        if "." in formula or formula == "": continue
        found = False
        for el in excl_els:
            if el in formula: 
                print formula
                found = True
        if found: continue

        if [formula, groupnum] not in existset:
            existset.append([formula, groupnum])
            icsdat[icsd[5:11]] = [formula, groupnum]

print "length", len(icsdat)
for key, value in icsdat.items():
    print key, value[0], value[1]


json.dump(icsdat, open("icsd.json", 'w'))
