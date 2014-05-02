import commands
import json
import pickle
import numpy as np
from collections import Counter, defaultdict
from pylab import *
from matplotlib.ticker import FuncFormatter
from ase.units import Bohr
from ase.utils import gcd
from atomic_constants import pauling, radius, Zval, Eatom, Emadelung, charge, mus
from sklearn import neighbors

Exptvol = json.load(open("exptvol.json",'r'))

class Atoms:
    def __init__(self, item=None):
        if item is not None:
            self.masses = np.array(item["atommasses_amu"])
            self.Z = np.array(item["atomvalences"])
            self.names = item["atomnames"]
#            self.Z = []
#            for name in self.names:
#                self.Z.append(Zval[name])
            self.positions = np.array(item["finalcartposmat"])
            self.natoms = int(item["numatom"])
            self.cell = np.array(item["finalbasismat"])
            self.formula = item["formula"]
            self.ncell = get_number_of_primitive_cell(item["atommasses_amu"])
            self.Eref = float(item["energyperatom"])
#            for name in self.names:
#                self.Eref -= mus[name] / self.natoms
#            self.eigenmat = np.array(item["eigenmat"])
            icsdstr = "{0:06d}".format(int(item["icsdnum"]))
            self.icsdno = icsdstr
#            name = self.names[np.argsort(self.masses)[0]]
#            self.spacegroup, self.exptvol = get_spacegroup_and_volume(name, icsdstr, self.natoms)
            self.calcvol = float(item["finalvolume_ang3"]) / self.natoms #self.ncell

def get_number_of_primitive_cell(Z):
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

def get_spacegroup_and_volume(name, icsdno, natoms):
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
#    M = set_all_coulumb_matrix(mset)
#    print "Finished coulomb matrix"

    # vladan's version
    pkl_file = open('Smatrix.pkl', 'rb')
    s_matrix = pickle.load(pkl_file)

    K = np.zeros((nset, nset))
    for i in range(nset):
        M1 = s_matrix["%s"%(mset[i].icsdno)]
        for j in range(nset):
            M2 = s_matrix["%s"%(mset[j].icsdno)]
            d = np.trace(np.dot(M1-M2,M1-M2))
#            K[i, j] = get_kernel(distance(M[i, :], M[j, :]), sigma, kernel=kernel)
            K[i, j] = get_kernel(d, sigma, kernel=kernel)
        K[i, i] += lamda
    print "Finished kernel"

    pkl_file.close()

    alpha = np.dot(np.linalg.inv(K), Eref) # not sure about the order
#    return M, alpha
    return alpha


def get_kernel(d, sigma, kernel="gaussian"):
    if kernel == "gaussian":
        return np.exp(- d**2 / (2.*sigma**2))
    elif kernel == "laplacian":
        return np.exp(-np.abs(d)/sigma)
    else:
        print "kernel not defined"
        XX

def estimation(mtrain, Etrain, M, alpha, sigma, mcross=None, Ecross=None, kernel="laplacian"):
    pkl_file = open('Smatrix.pkl', 'rb')
    s_matrix = pickle.load(pkl_file)

    nj = len(mtrain)
    if mcross is not None:
        ni = len(mcross)
        Eref = Ecross
#        Mref = set_all_coulumb_matrix(mcross)
        mset = mcross
    else:
        ni = nj
        Eref = Etrain
#        Mref = M
        mset = mtrain
    MAE = 0
    for i in range(ni):
        Eest = 0 # estimation for set number i
        Mref = s_matrix["%s"%(mset[i].icsdno)]
        for j in range(nj):
            M = s_matrix["%s"%(mtrain[j].icsdno)]
            d = np.trace(np.dot(Mref-M,Mref-M))
            Eest += alpha[j] * get_kernel(d, sigma, kernel=kernel)
#            Eest += alpha[j] * get_kernel(distance(Mref[i, :], M[j, :]), sigma, kernel=kernel)
        MAE += np.abs(Eest - Eref[i])
#        print mset[i].formula, mset[i].natoms, mset[i].ncell, Eest, Eref[i], Eest - Eref[i]
    pkl_file.close()
    return MAE


def knn_regression(mtrain, mcross, n_ngh, ndim,scaling):
    Etrain = np.array(get_Eref(mtrain))
    Ecross = np.array(get_Eref(mcross))

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
    plot(model.predict(Xcross),Ecross, 'ok')
    show()
    return np.nansum(np.abs(model.predict(Xcross) - Ecross)) / len(Ecross) # MAE


def read_json(filename = "data_RS.json"):
    d = json.load(open(filename, 'r'))
    formulas = []
    mset = []
    for i, item in enumerate(d):
        atoms = Atoms(item)
        try:
            atoms.exptvol = Exptvol[atoms.icsdno]
        except:
            continue
        if atoms.formula not in formulas: 
            formulas.append(atoms.formula)
        else:
            continue
#        if len(mset) > 0 and atoms.formula == mset[-1].formula: continue
        mset.append(atoms)
    print "Size of dataset : ", len(mset)//5*5
    del d

    Eref = get_Eref(mset)
    index = np.argsort(Eref)
    msetnew = []
    for i in index:
        msetnew.append(mset[i])

    return msetnew[:len(msetnew)//5*5] # return set that is divisable by 5, since its 5 fold 

def get_Eref(mset):
    Eref = []
    for atom in mset:
        Eref.append(atom.Eref)
    return Eref


def choose_lamda_sigma(mtrain, mcross):
    Etrain = get_Eref(mtrain)
    Ecross = get_Eref(mcross)

    for sigma in (10,50): #np.linspace(1,5,4):
        for lamda in (0.01, 0.001, ): #np.linspace(0.5, 2.5, 4):
            alpha = regression(mtrain, Etrain, sigma=sigma, lamda=lamda)
            M = 0
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
    

def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] == True:
        return s + r'$\%$'
    else:
        return s + '%'


def plot_Eref(mset):
    Eref = get_Eref(mset)
    plot_all(Eref, "Total energy per atom (eV)")


def plot_all(d, xl):
    hist(d, 50, normed=False)
    xlabel(xl, fontsize=14)
#    formatter = FuncFormatter(to_percent)
#    gca().yaxis.set_major_formatter(formatter)
    show()

def plot_natoms(mset):
    natoms = []
    for atom in mset:
        natoms.append(atom.natoms)
    plot_all(natoms, "Number of atoms in simulated cell")
    

def plot_error_in_volume(mset):
    volerror = []
    for atoms in mset:
        if atoms.exptvol is not None:
            volerror.append((atoms.calcvol - atoms.exptvol) / atoms.exptvol * 100)
            if np.abs(volerror[-1]) > 50:
                print atoms.formula, atoms.icsdno, atoms.exptvol, atoms.calcvol, volerror[-1]

    plot_all(volerror,"Error in volume (%) : (Vcalc - Vexpt) / Vexpt * 100")

def get_unique_spacegroups(mset):
    spacegroups = set()
    for atoms in mset:
        spacegroups.add(atoms.spacegroup)
    return spacegroups

def get_unique_elements(mset):
    elements = defaultdict(int)
    for atoms in mset:
        tmp = set()
        for name in atoms.names:
            tmp.add(name)
        for name in tmp:
            elements[name] += 1
    return elements

def plot_elements(mset):
    a = get_unique_elements(mset)
    X = np.arange(len(a))
    bar(X, a.values(), align="center", width=0.5)
    xticks(X, a.keys())
    ymax = max(a.values()) + 20
    ylim(0, ymax)
    show()

def write_csv(strtype='general'):
    if strtype == "general":
        mset = read_json(filename = "data.json")
        f = open('general.csv', 'w')
        print >> f, "(formula, std.mass, sum.mass, std.elecneg, sum.elecneg, calcvol, std.radius, sum.radius, Ecoh"
        for atoms in mset:
            elecneg = []
            rad = []
            for name in atoms.names:
                elecneg.append(pauling[name])
                rad.append(radius[name])
            sum_elecneg = np.sum(elecneg) / len(elecneg)
            std_elecneg = np.std(elecneg) 
            sum_mass = np.sum(atoms.masses) / len(atoms.masses)
            std_mass = np.std(atoms.masses)
            sum_radius = np.sum(rad) / len(rad)
            std_radius = np.std(rad)
            print >>f, "%s, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f"%(atoms.formula, std_mass, sum_mass, std_elecneg, sum_elecneg, 
                                                                                     atoms.calcvol, std_radius, sum_radius, atoms.Eref)
    elif strtype == "RS":
        mset = read_json(filename = "data_RS.json")
        f = open('RS.csv', 'w')
        print >> f, "(name0, name1, el, sum.mass, dmass, sum.elecneg, delecneg, calcvol, sum.radius, draius, dpos, Emadelung, Ecoh"
        for atoms in mset:
            elecneg1 = pauling[atoms.names[0]]
            elecneg2 = pauling[atoms.names[1]]
            volscaled = atoms.calcvol**(1./3.) / (radius[atoms.names[0]] * radius[atoms.names[1]]) * 10**4
            delecneg = np.abs(elecneg1-elecneg2)
            sqrtneg = np.std([elecneg1, elecneg2]) #np.sqrt(elecneg1*elecneg2)
    #        Elatt = Emadelung["%s"%(atoms.icsdno)]
            for el in atoms.names:
                if el in charge.keys():
                    Elatt = Emadelung["%s"%(atoms.icsdno)] / charge[el]**2
                    break
            
            dmass = np.abs(atoms.masses[0] - atoms.masses[1])
            d = atoms.positions[0] - atoms.positions[1]
            dpos = np.sqrt(np.inner(d, d))
    
            print >>f, "%s, %s, %s, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f"%(atoms.names[0], atoms.names[1], el, atoms.masses[0] + atoms.masses[1], np.abs(atoms.masses[0] - atoms.masses[1]), (elecneg1 + elecneg2), delecneg, atoms.calcvol, radius[atoms.names[0]] + radius[atoms.names[1]], np.abs(radius[atoms.names[0]] - radius[atoms.names[1]]), dpos, Elatt, atoms.Eref)




if __name__ == "__main__":
    mset = read_json("data.json")
    mtest, mset = get_testset(mset)
    mtrain, mcross, mset = get_train_validation_set(mset)
    choose_lamda_sigma(mtrain, mcross)
#    for n in range(5,10):
#        for s1 in (0.01, 0.05, 0.1, 0.3):
#            for s2 in (0.005, 0.01, 0.02, 0.05):
#            print n, s1, knn_regression(mtrain, mcross, n, 3, [1, 0.01, 1])

#    for names in mset:
#        print names.formula, names.icsdno
#    plot_error_in_volume(mset)
#    plot_Eref(mset)
#    plot_natoms(mset)
#    print get_unique_spacegroups(mset)
#    plot_elements(mset)
#    print get_unique_elements(mset)
#    write_csv("RS")

    

