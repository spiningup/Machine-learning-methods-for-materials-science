import numpy as np
import json, pickle
import commands
from collections import Counter, defaultdict
from ase.utils import gcd
from atomic_constants import mus, Eatom

Exptvol = pickle.load(open("exptvol.pkl",'r'))
avg_cord = json.load(open("avg_cord.json", 'r'))


class Atoms:
    def __init__(self, item=None, getadd=False):
        if item is not None:
            self.masses = np.array(item["atommasses_amu"])
            self.Z = np.array(item["atomvalences"])
            self.names = item["atomnames"]
#            self.Z = []
#            for name in self.names:
#                self.Z.append(1)#Zval[name])
            self.positions = np.array(item["finalcartposmat"])
            self.natoms = int(item["numatom"])
            self.cell = np.array(item["finalbasismat"])
            self.formula = item["formula"]
            self.ncell = self.get_number_of_primitive_cell(item["atommasses_amu"])
            self.Eref = float(item["energyperatom"])
            for name in self.names:
                self.Eref -= Eatom[name] / self.natoms
#                self.Eref -= mus[name] / self.natoms
#            self.eigenmat = np.array(item["eigenmat"])
            icsdstr = "{0:06d}".format(int(item["icsdnum"]))
            self.icsdno = icsdstr
            self.exptvol = Exptvol[self.icsdno][6]
            self.latt_a, self.latt_b, self.latt_c = np.sort(Exptvol[self.icsdno][0:3])
            self.alpha, self.beta, self.gamma = Exptvol[self.icsdno][3:6]
            self.avg_cord = avg_cord[self.icsdno]

            if getadd:
                # get stuff not in json file
                # get spacegroup by its name and icsdno
                name = self.names[np.argsort(self.masses)[0]]
                self.spacegroup, self.exptvol = self.get_spacegroup_and_volume(name, icsdstr, self.natoms)

            # calculated volume per atom
            self.calcvol = float(item["finalvolume_ang3"]) / self.natoms #self.ncell


    def get_number_of_primitive_cell(self, Z):
        # input a list of atomic masses, and output number of primitive cells
        b = Counter(Z)
        nlist = [i[1] for i in b.items()] # number of atoms for each atom species
        if len(nlist) == 1:
            return nlist[0]
        else:
            ncell = gcd(nlist[0], nlist[1])
            if len(nlist) > 2:
                for i in range(2, len(nlist)):
                    ncell = gcd(ncell, nlist[i])
            return ncell

    
    def get_spacegroup_and_volume(self, name, icsdno, natoms):
        # get experimental volume per formula unit
        output =  commands.getoutput("grep 'space' /home/jyan/icsd/%s/icsd_%s.cif"%(name, icsdno))
        if "No such file or directory" in output:
            return None, None
        
        volume =  float(commands.getoutput("grep '_cell_volume' /home/jyan/icsd/%s/icsd_%s.cif"%(name, icsdno)).split()[1])
        formulaunit = float(commands.getoutput("grep '_cell_formula_units_Z' /home/jyan/icsd/%s/icsd_%s.cif"%(name, icsdno)).split() [1])
        chemformula = commands.getoutput("grep '_chemical_formula_sum' /home/jyan/icsd/%s/icsd_%s.cif"%(name, icsdno)).split()
        na = []
        for s in chemformula:
            x = 0
            for i in s:
                if i.isdigit():
                    x = int(i) + 10 * x
            na.append(x)
        if len(na) == 1:
            formulaunit *= na[0]
        else:
            ncell = gcd(na[0], na[1])
            if len(na) > 2:
                for i in range(2, len(na)):
                    ncell = gcd(ncell, na[i])
            formulaunit *= ncell
    
    #    try:
    #        volume = read.icsd_cif_a("/home/jyan/icsd/%s/icsd_%s.cif"%(name, icsdno)).volume.magnitude
    #    except:
    #        return None, None
    
        return output.split("'")[1], volume / formulaunit


if __name__ == "__main__":
    d = json.load(open("data_RS.json", 'r'))
    for item in d:
        atoms = Atoms(item)
        print atoms.icsdno, atoms.formula, atoms.spacegroup, atoms.calcvol * atoms.natoms/atoms.ncell, atoms.exptvol

    
