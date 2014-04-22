import numpy as np
from ase.units import Bohr
#from ase.io import read
import json
import random
from pylab import *

class Atoms:
    def __init__(self, item=None):
        if item is not None:
#            self.Z = np.array(item["atommasses_amu"])
            self.Z = np.array(item["atomvalences"])
            self.positions = np.array(item["finalcartposmat"])
            self.natoms = int(item["numatom"])
            self.cell = np.array(item["finalbasismat"])
            self.formula = item["formula"]
            self.Eref = float(item["energynoentrp"]) / self.natoms


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

def set_all_coulumb_matrix(mset):

    nset = len(mset) # number of training set
    max_natoms = 0   # maximum number of atoms in the training set
    for atoms in mset:
        if atoms.natoms > max_natoms:
            max_natoms = atoms.natoms
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
    else:
        ni = nj
        Eref = Etrain
        Mref = M
    MAE = 0
    for i in range(ni):
        Eest = 0 # estimation for set number i
        for j in range(nj):
            Eest += alpha[j] * get_kernel(distance(Mref[i, :], M[j, :]), sigma, kernel=kernel)
        MAE += np.abs(Eest - Eref[i])
#        print Eest, Eref[i], Eest - Eref[i]
    return MAE

def read_json(filename = "data.json"):
    d = json.load(open(filename, 'r'))
    mset = []
    for i, item in enumerate(d):
        atoms = Atoms(item)
        if i > 0 and atoms.formula == mset[-1].formula: continue
        mset.append(atoms)
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

    for sigma in (15, 20): #np.linspace(1,5,4):
        for lamda in (0.1, ): #np.linspace(0.5, 2.5, 4):
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
#    Eref = get_Eref(mset)
#    for i in range(len(Eref)):
#        if Eref[i] < -10:
#            print mset[i].formula, Eref[i]
#    from pylab import *
#    hist(Eref, 20)
#    show()
