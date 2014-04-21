import numpy as np
from ase.units import Bohr
#from ase.io import read
import json
import random

class Atoms:
    def __init__(self, item=None):
        if item is not None:
            self.Z = np.array(item["atommasses_amu"])
            self.positions = np.array(item["finalcartposmat"])
            self.natoms = int(item["numatom"])
            self.cell = np.array(item["finalbasismat"])
            self.formula = item["formula"]


def set_coulumb_matrix(atoms):
    
    na = atoms.natoms
    Z = atoms.Z #get_atomic_numbers()
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
    
    return E

def distance(M1, M2):
    n = max(len(M1), len(M2))
    M1 = np.append(M1, np.zeros(n-len(M1)))
    M2 = np.append(M2, np.zeros(n-len(M2)))
    return np.sqrt(np.dot(M1-M2, M1-M2))

def regression(mset, Eref, sigma, lamda):
    # lamda for regularization, sigma for gaussian damping
    nset = len(mset) # number of training set
    max_natoms = 0   # maximum number of atoms in the training set
    for atoms in mset:
        if atoms.natoms > max_natoms:
            max_natoms = atoms.natoms
    M = np.zeros((nset, max_natoms))
    for i, atoms in enumerate(mset):
        Mtmp = set_coulumb_matrix(atoms)
        M[i,:len(Mtmp)] = Mtmp[:]
    print "Finished coulomb matrix"

    K = np.zeros((nset, nset))
    for i in range(nset):
#        print i
        for j in range(nset):
            K[i, j] = np.exp(- distance(M[i, :], M[j, :])**2 / (2. * sigma**2))
        K[i, i] += lamda
    print "Finished kernel"

    alpha = np.dot(np.linalg.inv(K), Eref) # not sure about the order
    return M, alpha

def estimation(mset, Eref, alpha, M, sigma):
    nset = len(mset)
    MAE = 0
    for i in range(nset):
        Eest = 0 # estimation for set number i
        for j in range(nset):
            Eest += alpha[j] * np.exp(-distance(M[i, :], M[j, :])**2 / (2. * sigma**2)) 
        MAE += np.abs(Eest - Eref[i])
#        print mset[i].formula, Eest, Eref[i], Eest - Eref[i]
    return MAE

def read_json(filename = "data.json"):
    d = json.load(open(filename, 'r'))
    Eref = []
    mset = []
    for i, item in enumerate(d):
        atoms = Atoms(item)
        if i > 0 and atoms.formula == mset[-1].formula: continue
        mset.append(atoms)
        Eref.append(float(item["energyperatom"])) #* atoms.natoms
    del d

    return mset, Eref

def choose_lamda_sigma():
    mset, Eref = read_json()
    for sigma in (0.1, ): #np.linspace(1,5,4):
        for lamda in (0.005,): #np.linspace(0.5, 2.5, 4):
            M, alpha = regression(mset, Eref, sigma=sigma, lamda=lamda)
            MAE = estimation(mset, Eref, sigma=sigma, alpha=alpha, M=M)
            print sigma, lamda, min(alpha), max(alpha), MAE / len(mset)


if __name__ == "__main__":
    mset, Eref = read_json()
    stratified_sampling(mset, Eref)

# examine histgram
    for i in range(len(Eref)):
        if Eref[i] < -10:
            print mset[i].formula, Eref[i]
    from pylab import *
    hist(Eref, 20)
    show()
