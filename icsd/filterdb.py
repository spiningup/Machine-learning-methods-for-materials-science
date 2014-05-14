import numpy as np
import json
import os
from collections import defaultdict
from ase.utils import gcd

def filter_elements(d, el="O", incl=True):
    dnew = {}
    for key, item in d.items():
        if el in item[0] and incl:
            dnew[key] = item
        if el not in item[0] and (not incl):
            dnew[key] = item
    return dnew

def filter_nel(d, l="2"): # filter number of elements
    dnew = {}
    for key, item in d.items():
        formula = item[0].split()
        if len(formula) == l:
            dnew[key] = item
    return dnew

def filter_calculated(d, d0):
    dnew = {}
    for key, item in d.items():
        if key in d0:
            print key
            continue
        else:
            dnew[key] = item
    return dnew

def filter_disorder(d):
    dnew = {}
    for key, item in d.items():
        formula = item[0].replace("'",'').split()
        n = []
        for el in formula:
            for i in el:
                if not i.isdigit():
                    el = el.replace(i, '')
            n.append(int(el))

        ncell = n[0]
        if len(n) >= 2:
            ncell = gcd(n[0], n[1])
            if len(n) > 2:
                for i in range(2, len(n)):
                    ncell = gcd(ncell, n[i])

        natoms = sum(n) / ncell
        if item[2] != 0 and (item[2] % natoms == 0 or natoms % item[2] == 0):
            dnew[key] = item
        else:
            print key, formula, item[2], natoms
            
    return dnew
                
def filter_natoms(d, nmin=0, nmax=100):
    dnew = {}
    for key, item in d.items():
        if item[2] <= nmax and item[2] > nmin:
            dnew[key] = item
    return dnew

def formula_with_multiple_spacegroups(d):
    dnew = defaultdict(list)
    for key, item in d.items():
        formula, groupnum, natoms = item[0], item[1], item[2]
        dnew[formula].append([key, groupnum, natoms]) # icsdno, spacegroup, natoms

    dnew2 = {key: value for key, value in dnew.items() if len(value) > 1}
    return dnew2

def add_otherspacegroups(d, dnew):
    for key, item in d.items():
        formula = item[0]
        if formula in dnew.keys(): # if the formula has multiple spacegroups, check whether other spacegroups are in the dataset
            for i in dnew[formula]:
                if i[0] not in d.keys(): # icsd not in dataset
                    d[i[0]] = [formula, i[1], i[2]] # add spacegroup and natoms
    return d

def filter_magnetic_elements(d):
    mag = ['V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 
           'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag']

    dmag = {}; dnonmag = {}
    for key, item in d.items():
        nmag = 0
        for el in mag:
            if el in item[0]:
                nmag += 1
                dmag[key] = item
        if key not in dmag.keys():
            dnonmag[key] = item
        else:
            dmag[key] = item + [nmag,]
    return dnonmag, dmag

if __name__ == "__main__":

    d = json.load(open('icsd_wi_natoms.json', 'r'))
    dtmp = json.load(open('icsdnum_NRELDB.json', 'r'))
    d0 = []
    for i in dtmp:
        icsdno = i.values()[0]
        icsdstr = "{0:06d}".format(int(icsdno))
        d0.append(icsdstr)

#    d = filter_disorder(d)
#    d = filter_calculated(d, d0)
#    json.dump(d, open("icsd.json", 'w'))
    
    d = json.load(open('icsd.json', 'r'))
    for i in ([0, 10, 20, 30, 40, 50, 60, 70, 80, 90]):
        print i, "to", i+10, len(filter_natoms(d, nmin=i, nmax=i+10))

#    d = filter_elements(d,"O")
#    for key, value in d.items():
#        print key, value


