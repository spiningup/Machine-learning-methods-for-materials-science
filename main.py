import numpy as np
from ase.units import Bohr
import json
import random
from pylab import *
from ase.utils import gcd
from collections import Counter

Zval = {
"Al" : 3, 
"As" : 3, #8-3,
"Ba" : 2, 
"Bi" : 3, #8-3,
"Ca" : 2, 
"Cl" : 1, #8-1, 
"F"  : 1, #8-1, 
"Ga" : 3, 
"Ge" : 4, 
"In" : 3, 
"K"  : 1,  
"Mg" : 2, 
"N"  : 3, #8-3, 
"Na" : 1, 
"O"  : 2, #8-2,
"P"  : 3, #8-3, 
"Rb" : 1, 
"S"  : 2, #8-2, 
"Sb" : 3, #8-3, 
"Se" : 2, #8-2,
"Si" : 4,
"Sn" : 4, 
"Sr" : 2, 
"Te" : 2, #8-2,
"Y"  : 3,
}

class Atoms:
    def __init__(self, item=None):
        if item is not None:
#            self.Z = np.array(item["atommasses_amu"])
            self.Z = np.array(item["atomvalences"])
            self.names = item["atomnames"]
            self.Z = []
            for name in self.names:
                self.Z.append(Zval[name])
            self.positions = np.array(item["finalcartposmat"])
            self.natoms = int(item["numatom"])
            self.cell = np.array(item["finalbasismat"])
            self.formula = item["formula"]
            self.ncell = get_number_of_primitive_cell(item["atommasses_amu"])
            self.Eref = float(item["energynoentrp"]) / self.natoms
#            self.eigenmat = np.array(item["eigenmat"])

def get_number_of_primitive_cell(Z):
    b = Counter(Z)
    nlist = [i[1] for i in b.items()] # number of atoms for each atom species
    if len(nlist) == 1:
        return nlist[0]
    else:
        ncell = gcd(nlist[0], nlist[1])
        if len(nlist) > 2:
            for i in range(2, len(nlist)):
                ncell = gcd(ncell, nlist[2])
        return ncell

def set_coulumb_matrix(atoms):
    
    na = atoms.natoms
    Z = atoms.Z 
    V = np.zeros((na, na))
#    n1, n2, n3 = 2,2,2    

    for i in range(na):
        for j in range(na):
            if i == j:
                V[i, j] = 0.5 * Z[i]**2.4
            else:
                d = (atoms.positions[i] - atoms.positions[j]) / Bohr
                V[i, j] = Z[i] * Z[j] / np.sqrt(np.dot(d, d))

# ------- neighbors -------------
#            for i1 in range(-n1/2, n1/2+1):
#                for i2 in range(-n2/2, n2/2+1):
#                    for i3 in range(-n3/2, n3/2+1):
#                        if i1 == i2 == i3 == 0:
#                            continue
#                        rb = atoms.positions[j] + np.dot(np.array([i1, i2, i3]), atoms.cell)
#                        d = (atoms.positions[i] - rb) / Bohr
#                        V[i, j] += Z[i] * Z[j] / np.sqrt(np.dot(d, d)) #* np.exp(-np.dot(d, d)/10)

    E = np.linalg.eig(V)[0] / na**2

    #E = atoms.eigenmat[0][0]
    
    return E

def set_all_coulumb_matrix(mset):

    nset = len(mset) # number of training set
    max_natoms = 0   # maximum number of atoms in the training set
    for atoms in mset:
        if atoms.natoms > max_natoms:
            max_natoms = atoms.natoms
#    for atoms in mset:
#        if len(atoms.eigenmat[0][0]) > max_natoms:
#            max_natoms = len(atoms.eigenmat[0][0])
    M = np.zeros((nset, max_natoms))
    for i, atoms in enumerate(mset):
        Mtmp = set_coulumb_matrix(atoms)
        M[i,:len(Mtmp)] = Mtmp[:]
    return M



def distance(M1, M2):
    n = max(len(M1), len(M2))
    M1 = np.append(M1, np.zeros(n-len(M1)))
    M2 = np.append(M2, np.zeros(n-len(M2)))
    return np.sqrt(np.dot(M1-M2, M1-M2))

def regression(mset, Eref, sigma, lamda, kernel="laplacian"):
    # lamda for regularization, sigma for gaussian damping
    nset = len(mset) # number of training set
    M = set_all_coulumb_matrix(mset)
    print "Finished coulomb matrix"

    K = np.zeros((nset, nset))
    for i in range(nset):
        for j in range(nset):
            K[i, j] = get_kernel(distance(M[i, :], M[j, :]), sigma, kernel=kernel)
        K[i, i] += lamda
    print "Finished kernel"

    alpha = np.dot(np.linalg.inv(K), Eref) # not sure about the order
    return M, alpha


def get_kernel(d, sigma, kernel="gaussian"):
    if kernel == "gaussian":
        return np.exp(- d**2 / (2.*sigma**2))
    elif kernel == "laplacian":
        return np.exp(-np.abs(d)/sigma)
    else:
        print "kernel not defined"
        XX

def estimation(mtrain, Etrain, M, alpha, sigma, mcross=None, Ecross=None, kernel="laplacian"):
    nj = len(mtrain)
    if mcross is not None:
        ni = len(mcross)
        Eref = Ecross
        Mref = set_all_coulumb_matrix(mcross)
        mset = mcross
    else:
        ni = nj
        Eref = Etrain
        Mref = M
        mset = mtrain
    MAE = 0
    for i in range(ni):
        Eest = 0 # estimation for set number i
        for j in range(nj):
            Eest += alpha[j] * get_kernel(distance(Mref[i, :], M[j, :]), sigma, kernel=kernel)
        MAE += np.abs(Eest - Eref[i])
#        print mset[i].formula, mset[i].natoms, mset[i].ncell, Eest, Eref[i], Eest - Eref[i]
    return MAE

def read_json(filename = "data.json"):
    d = json.load(open(filename, 'r'))
    mset = []
    for i, item in enumerate(d):
        atoms = Atoms(item)
        if len(mset) > 0 and atoms.formula == mset[-1].formula: continue
        mset.append(atoms)
    print "Size of dataset : ", len(mset)//5*5
    del d

    return mset[:len(mset)//5*5] # return set that is divisable by 5, since its 5 fold 

def get_Eref(mset):
    Eref = []
    for atom in mset:
        Eref.append(atom.Eref)
    return Eref

def plot_histgram(mset):
    Eref = get_Eref(mset)
    hist(Eref, 20)
    show()


def choose_lamda_sigma(mtrain, mcross):
    Etrain = get_Eref(mtrain)
    Ecross = get_Eref(mcross)

    for sigma in (0.1, ): #np.linspace(1,5,4):
        for lamda in (0.01, 0.001, 0.0001): #np.linspace(0.5, 2.5, 4):
            M, alpha = regression(mtrain, Etrain, sigma=sigma, lamda=lamda)
            MAEtrain =  estimation(mtrain, Etrain, M, alpha, sigma)
            MAEcross = estimation(mtrain, Etrain, M, alpha, sigma, mcross, Ecross)

            print sigma, lamda, MAEtrain / len(mtrain), MAEcross / len(mcross)


def get_testset(mset):
    mtest = []
    i = 0
    while i < len(mset):
        idx = np.random.randint(0, 5)
        mtest.append(mset[i+idx])
        mset.pop(i+idx)
        i += 4
        
    return mtest, mset

def get_train_validation_set(mset):
    mcross = []
    mtrain = mset[:]
    i = 0
    while i < len(mtrain):
        idx = np.random.randint(0, 4)
        mcross.append(mtrain[i+idx])
        mtrain.pop(i+idx)
        i += 3
        
    return mtrain, mcross, mset # mset is a sum of mtrain and mcross
    
    

if __name__ == "__main__":
    mset = read_json()
    mtest, mset = get_testset(mset)
    mtrain, mcross, mset = get_train_validation_set(mset)
    choose_lamda_sigma(mtrain, mcross)

# examine coulumb matrix
#    M = set_all_coulumb_matrix(mset)
#    for i in range(len(M)):
#        print i
#        print M[i, :]
    
# examine histgram
#    natoms_all = []
#    for atoms in mset:
#        natoms_all.append(atoms.natoms)
#    hist(natoms_all, 20)
#    show()
#    Eref = get_Eref(mset)
#    for i in range(len(Eref)):
#        if Eref[i] < -10:
#            print mset[i].formula, Eref[i]
#    from pylab import *
#    hist(Eref, 20)
#    show()
