import numpy as np
from read_json import read_json, get_elements_map
from split_dataset import *
from ml import *
from pylab import *
np.random.seed(0)

mset = read_json("larger_data/data.json", energytype="formation")
mcross, mtrain = get_testset(mset)
mcross = read_json("CoSi.json", energytype="formation")
elmap = get_elements_map(mset)
Xtrain, Xcross, Etrain, Ecross = get_X(mtrain, mcross, elmap=elmap, elmethod="composition")
    
alpha = pickle.load(open('alpha.pkl', 'r'))

sigma = 12
print mcross[0].formula
for i in range(len(Ecross)):
    Eest = 0 # estimation for set number i
    Econ = []
    ECo = 0
    for j in range(len(Etrain)):
        d = Xtrain[j] - Xcross[i]
        dd = np.sqrt(np.inner(d, d))
        Etmp = alpha[j] * get_kernel(dd, sigma, kernel="gaussian")
        Eest += Etmp
        Econ.append(np.abs(Etmp))
        if "Co" in mtrain[j].formula:
            ECo += Etmp
    print Eest, ECo

    XX
    idx = np.argsort(Econ)
    for k in range(len(Etrain)):
        if "Co" in mtrain[idx[k]].formula:
            print mtrain[idx[k]].formula, Econ[idx[k]]
    XX
    print Eest
