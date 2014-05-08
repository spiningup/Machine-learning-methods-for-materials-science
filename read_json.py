import json
import numpy as np
from collections import defaultdict
from setup_atoms import Atoms, read_helper_data


def read_json(filename = "data.json", energytype="atomization"):
    d = json.load(open(filename, 'r'))
    formulas = []
    mset = []
    if filename.split("/")[0] == filename:
        prefix = "."
    else:
        prefix = filename.split("/")[0]
    Exptvol, cord, coulomb = read_helper_data(prefix)
    for i, item in enumerate(d):
        atoms = Atoms(item, Exptvol, cord, coulomb, energytype)
        if atoms.Eref is None: continue
#        if len(atoms.formula.split()) <=2 : continue
#        if "La" in atoms.names: continue
#        if "Y" in atoms.names: continue

        # ignore formula already in dataset
        formula = " ".join(sorted(atoms.formula.split()))
        atoms.formula = formula
        if formula not in formulas: 
            formulas.append(formula)
        else:
            continue

        volerror = np.abs(atoms.calcvol - atoms.exptvol) / atoms.exptvol 
        if volerror > 0.2: continue

        mset.append(atoms)
    del d

    # rank the dataset by Eref
    Eref = attribute_tolist(mset, attr="Eref")
    index = np.argsort(Eref)
    msetnew = []
    for i in index:
        msetnew.append(mset[i])

    print "Size of dataset : ", len(msetnew)//5*5

    return msetnew[:len(msetnew)//5*5] # return set that is divisable by 5, since its 5 fold 

def attribute_tolist(mset, attr="Eref", unique=False):
    if not unique:
        attrlist = []
    else:
        attrlist = set()
    for atom in mset:
        if not unique:
            attrlist.append(getattr(atom, attr))
        else:
            attrlist.add(getattr(atom, attr))
        
    return attrlist

def get_unique_elements(mset):
    elements = defaultdict(int)
    for atoms in mset:
        tmp = set()
        for name in atoms.names:
            tmp.add(name)
        for name in tmp:
            elements[name] += 1
    return elements

def get_elements_map(mset):
    elements = get_unique_elements(mset)
    elmap = {}
    for i, el in enumerate(elements):
        elmap[el] = i
    return elmap

if __name__ == "__main__":
    mset = read_json("data.json")
    Eref = attribute_tolist(mset, attr="Eref", unique=False)
    elements = get_unique_elements(mset)
    
#    print Eref
    print elements, len(elements)

    numel_performula = defaultdict(int)
    for atoms in mset:
        l = len(atoms.formula.split())
        numel_performula[l] += 1
    print numel_performula
#        print atoms.formula, atoms.Eref
