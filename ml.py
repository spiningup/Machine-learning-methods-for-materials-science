from sklearn import neighbors
import numpy as np
from ase.units import Bohr
from read_json import attribute_tolist, read_json
from split_dataset import *
from atomic_constants import pauling, radius, Zval, Eatom, Emadelung, charge, mus

class kernel_ridge_regression:
    def __init__(self, mtrain, mcross, lamda, sigma, matrixtype=1):
        self.choice = matrixtype
        self.sigma = sigma
        self.lamda = lamda
        self.mtrain = mtrain
        self.mcross = mcross
        

    def set_coulumb_matrix(self, atoms):
        # set coulomb matrix for a single solid
        na = atoms.natoms
        Z = atoms.Z 
        V = np.zeros((na, na))
    
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
    
    def set_all_coulumb_matrix(self, choice=1, mset=None):
    
        if choice == 1:
            nset = len(mset) # number of training set
            M = {} 
            for i, atoms in enumerate(mset):
                M["%s"%(atoms.icsdno)] = self.set_coulumb_matrix(atoms)
            return M # matrix dictionary by icsdno
    
        elif choice == 2:
            pkl_file = open('Smatrix.pkl', 'rb')
            s_matrix = pickle.load(pkl_file)
            pkl_file.close()
            return s_matrix
        else:
            "Not implemented"
            XX
    
    
    def distance(self, M1, M2):
        if M1.ndim == 1 and M2.ndim == 1:
            n = max(len(M1), len(M2))
            M1 = np.append(M1, np.zeros(n-len(M1)))
            M2 = np.append(M2, np.zeros(n-len(M2)))
            return np.sqrt(np.dot(M1-M2, M1-M2))
        elif M1.ndim == 2 and M2.ndim == 2:
            return np.trace(np.dot(M1-M2,M1-M2))
        else:
            XX

    
    def regression(self, mset, M, sigma, lamda, kernel="laplacian"):
        # lamda for regularization, sigma for gaussian damping
        nset = len(mset) # number of training set
        Eref = attribute_tolist(mset, attr="Eref")
    
        K = np.zeros((nset, nset))
        for i in range(nset):
            M1 = M["%s"%(mset[i].icsdno)]
            for j in range(nset):
                M2 = M["%s"%(mset[j].icsdno)]
                K[i, j] = self.get_kernel(self.distance(M1, M2), sigma, kernel=kernel)
            K[i, i] += lamda
        print "Finished kernel"
    
        alpha = np.dot(np.linalg.inv(K), Eref) # not sure about the order
        return alpha
    
    
    def get_kernel(self, d, sigma, kernel="gaussian"):
        if kernel == "gaussian":
            return np.exp(- d**2 / (2.*sigma**2))
        elif kernel == "laplacian":
            return np.exp(-np.abs(d)/sigma)
        else:
            print "kernel not defined"
            XX
    
    def estimation(self, mtrain, M, alpha, sigma, mcross=None, Mref=None, kernel="laplacian"):
    
        ni = len(mcross)
        nj = len(mtrain)
        Eref = attribute_tolist(mcross, attr="Eref")

        MAE = 0
        for i in range(ni):
            Eest = 0 # estimation for set number i
            M1 = Mref["%s"%(mcross[i].icsdno)]
            for j in range(nj):
                M2 = M["%s"%(mtrain[j].icsdno)]
                Eest += alpha[j] * self.get_kernel(self.distance(M1, M2), sigma, kernel=kernel)
            MAE += np.abs(Eest - Eref[i])
    #        print mset[i].formula, mset[i].natoms, mset[i].ncell, Eest, Eref[i], Eest - Eref[i]
        return MAE / ni

    def run(self, sigma=None, lamda=None):
        if sigma is None and lamda is None:
            sigma = self.sigma
            lamda = self.lamda
        M = self.set_all_coulumb_matrix(choice=self.choice, mset=self.mtrain)
        Mref = self.set_all_coulumb_matrix(choice=self.choice, mset=self.mcross)
        alpha = self.regression(self.mtrain, M, sigma, lamda)
        MAEtrain =  self.estimation(self.mtrain, M, alpha, sigma, self.mtrain, M)
        MAEcross = self.estimation(self.mtrain, M, alpha, sigma, self.mcross, Mref)
        
        return MAEtrain, MAEcross


    def choose_lamda_sigma(self, sigmalist, lamdalist):
        for sigma in sigmalist: 
            for lamda in lamdalist:
                print sigma, lamda, self.run(sigma, lamda)
    

def knn_regression(mtrain, mcross, n_ngh, ndim,scaling):
    Etrain = np.array(attribute_tolist(mtrain, attr="Eref"))
    Ecross = np.array(attribute_tolist(mcross, attr="Eref"))
    
    Xtrain = np.zeros((len(Etrain), ndim)); Xcross = np.zeros((len(Ecross), ndim))
    for mset in (mtrain, mcross):
        for i, atoms in enumerate(mset):
            elecneg = 0
            rad = 0
            for name in atoms.names:
                elecneg += pauling[name] 
                rad += radius[name]
            if mset == mtrain:
                Xtrain[i] = atoms.exptvol**(1./3.), rad/atoms.natoms, elecneg/atoms.natoms
            elif mset == mcross:
                Xcross[i] = atoms.exptvol**(1./3.), rad/atoms.natoms, elecneg/atoms.natoms
    for i in range(ndim):
        Xtrain[:,i] *= scaling[i]
        Xcross[:,i] *= scaling[i]
    
    n_neighbors = n_ngh
    knn = neighbors.KNeighborsRegressor(n_neighbors, weights="distance")
    model = knn.fit(Xtrain, Etrain)

    return np.nansum(np.abs(model.predict(Xcross) - Ecross)) / len(Ecross) # MAE


if __name__ == "__main__":
    mset = read_json("data.json")
    mtest, mset = get_testset(mset)
    mtrain, mcross, mset = get_train_validation_set(mset)
    sigma = 50
    lamda = 0.01

    kRR = kernel_ridge_regression(mtrain, mcross, lamda, sigma, matrixtype=1)
    print kRR.run()
    kRR.choose_lamda_sigma([10, 50], [0.01, 0.001])

    print knn_regression(mtrain, mcross, 5, 3, [1, 0.01, 1])

