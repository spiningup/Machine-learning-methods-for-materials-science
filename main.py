import pickle
import numpy as np
from atomic_constants import pauling, radius, Zval, Eatom, Emadelung, charge, mus
from read_json import read_json
from split_dataset import *
from ml import *

mset = read_json("data.json")
mtest, mset = get_testset(mset)
mtrain, mcross, mset = get_train_validation_set(mset)
sigma = 50
lamda = 0.01
    
#kRR = kernel_ridge_regression(mtrain, mcross, lamda, sigma, matrixtype=1)
#print kRR.run()

#kRR = kernel_ridge_regression(mtrain, mcross, lamda, sigma, matrixtype=2)
#kRR.choose_lamda_sigma([10, 50], [0.01, 0.001])

# test knn
#for kernel in ([None, 'rbf', 'poly','cosine']):
#    for scaling in ([1, 2, 3]):
#        for weights in (["distance", "uniform"]):
#            print kernel, scaling, weights, knn_regression(mtrain, mcross, 5, kernel=kernel, scaling=scaling, weights=weights)

for metric in (["euclidean", "manhattan", "chebyshev", "minkowski",]):
    print metric, knn_regression(mtrain, mcross, 5, metric=metric)
#print krr_regression(mtrain, mcross, 50, 0.01)

#pca_decomposition(mtrain, mcross)





    

