import numpy as np
import json
import sys
from read_json import read_json, get_elements_map
from split_dataset import *
from ml import *
from pylab import *
np.random.seed(0)

def write_json(el1, el2):
    data = []
    data.append({"formula": "%s 1"%(el1), "FERE": 0})
    for i in range(9):
        data.append({"formula": "%s %d %s %d"%(el1, 9-i, el2, i+1), "FERE": 0})
    data.append({"formula": "%s 1"%(el2), "FERE": 0})
    json.dump(data, open("%s%s.json"%(el1, el2), 'w'))

elmethod = "composition"
sigma = 12 ; lamda = 0.0001 ; kernel = "gaussian"

MAEtrain = []
MAEcross = []
el1 = sys.argv[1]
el2 = sys.argv[2]

maxrun = 5
avgEpredict = np.zeros(11)
for irun in range(maxrun):
    write_json(el1, el2)
#    mset = read_json("include_ML_natoms_30/data.json", energytype="formation")
    mset = read_json("tests/data.json", energytype="formation")
    mcross, mtrain = get_testset(mset)
    mcross = read_json("%s%s.json"%(el1, el2), energytype="formation")
    #mtest, mset = get_testset(mset)
    #mtrain, mcross, mset = get_train_validation_set(mset)
    elmap = get_elements_map(mset)

    result = krr_regression(mtrain, mcross, sigma, lamda, kernel=kernel, elmap=elmap, elmethod=elmethod,
                            loadalpha=False, alphanum=irun)
    MAEtrain.append(result[0])
    MAEcross.append(result[1])
    Epredict = result[2]
    print result[0], result[1]

    print "element formation energy :", Epredict[0], Epredict[-1]
#    f = open("%sSi_result.dat"%(el), 'a')
#    for i in range(len(Epredict)):
#        print >> f, Epredict[i]
#    f.close()

    print len(Epredict)
    for i in range(1, 10):
        Epredict[i] = Epredict[i] - Epredict[-1]*i*0.1 - Epredict[0]*(10-i)*0.1
    Epredict[0] = Epredict[-1] = 0

    avgEpredict += np.array(Epredict)
    plot(np.arange(11)*0.1, avgEpredict/(irun+1), '-ok')
    show()
