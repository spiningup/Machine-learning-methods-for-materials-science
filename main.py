import numpy as np
from read_json import read_json, get_elements_map
from split_dataset import *
from ml import *
from pylab import *
#np.random.seed(0)


elmethod = "composition"
sigma = 12 ; lamda = 0.0001 ; kernel = "gaussian"
maxrun = 20

#for sigma in (10, 11, 12, 13, 14):
#    print krr_regression(mtrain, mcross, sigma, lamda, kernel=kernel, elmap=elmap, elmethod=elmethod)

MAEtrain = []
MAEcross = []
for i in range(maxrun):
    mset = read_json("data.json")
    mcross, mtrain = get_testset(mset)
    #mtest, mset = get_testset(mset)
    #mtrain, mcross, mset = get_train_validation_set(mset)
    elmap = get_elements_map(mset)

    result = krr_regression(mtrain, mcross, sigma, lamda, kernel=kernel, elmap=elmap, elmethod=elmethod)
    MAEtrain.append(result[0])
    MAEcross.append(result[1])
    print result

hist(MAEtrain, 5)
hist(MAEcross, 5)
show()


