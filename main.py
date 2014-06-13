import numpy as np
from read_json import *
from split_dataset import *
from ml import *
from visualize import *
from pylab import *
#np.random.seed(0)


elmethod = "composition"
sigma = 13 ; lamda = 0.001 ; kernel = "gaussian"
maxrun = 1

MAEtrain = []
MAEcross = []
for i in range(maxrun):
    mset = read_json("include_ML_natoms_30/data.json", energytype="formation")
    mcross, mtrain = get_testset(mset)
#    mcross = read_json("MnSi.json", energytype="formation")
    #mtest, mset = get_testset(mset)
    #mtrain, mcross, mset = get_train_validation_set(mset)
    elmap = get_elements_map(mset)

    result = krr_regression(mtrain, mcross, sigma, lamda, kernel=kernel, elmap=elmap, elmethod=elmethod,
                            loadalpha=False, alphanum=i)
    MAEtrain.append(result[0])
    MAEcross.append(result[1])
    print result[0], result[1]

#hist(MAEtrain, 10)
#hist(MAEcross, 10)
#show()

for i in range(5, 6):
    mset = read_json("include_ML_natoms_30/data.json", energytype="formation")
#    mset = read_json("larger_data/data.json", energytype="atomization")
           
#    msetnew = []
#    for atoms in mset:
#        if "O" in atoms.formula: continue# and j == 0: pass
#        if len(atoms.formula.split()) > 2: continue
#        else: msetnew.append(atoms)
#    mset = msetnew
#    print len(mset)
#    for atoms in mset:
#        print atoms.formula
#    plot_dict(get_unique_elements(mset))
    mset = get_partial_dataset(mset, i)
    print i, len(mset)
    elmap = get_elements_map(mset)
    mcross, mtrain = get_testset(mset)
#mtest, mset = get_testset(mset)
#mtrain, mcross, mset = get_train_validation_set(mset)

    print "knn MAE", knn_regression(mtrain, mcross, 5, selectf=False) 
    print "knn MAE", knn_regression(mtrain, mcross, 5, elmap=elmap, elmethod="composition") 
#    print "knn MAE", knn_regression(mtrain, mcross, 5, elmap=elmap, elmethod="coordination") 
#    print "knn MAE", knn_regression(mtrain, mcross, 5, elmap=elmap, elmethod="inverse_cord") 
#    print "knn MAE", knn_regression(mtrain, mcross, 5, elmap=elmap, elmethod="coulomb_ZiZj/d") 
