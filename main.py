import pickle
import numpy as np
from atomic_constants import pauling, radius, Zval, Eatom, Emadelung, charge, mus
from read_json import read_json, get_elements_map
from split_dataset import *
from ml import *

mset = read_json("data.json")
mtest, mset = get_testset(mset)
mtrain, mcross, mset = get_train_validation_set(mset)
elmap = get_elements_map(mset)

#sigma = 50
#lamda = 0.01    
#kRR = kernel_ridge_regression(mtrain, mcross, lamda, sigma, matrixtype=1)
#print kRR.run()

#kRR = kernel_ridge_regression(mtrain, mcross, lamda, sigma, matrixtype=2)
#kRR.choose_lamda_sigma([10, 50], [0.01, 0.001])

# test knn, kernel=None, scaling=1, weights=distance, metric="minkowsk" are the best options
#for kernel in ([None, 'rbf', 'poly','cosine']):
#    for scaling in ([1, 2, 3]):
#        for weights in (["distance", "uniform"]):
#            print kernel, scaling, weights, knn_regression(mtrain, mcross, 5, kernel=kernel, scaling=scaling, weights=weights)
#for metric in (["euclidean", "manhattan", "chebyshev", "minkowski",]):
#    print metric, knn_regression(mtrain, mcross, 5, metric=metric)
# select features
#print "knn MAE", knn_regression(mtrain, mcross, 5, selectf=True) 

#print "knn MAE", knn_regression(mtrain, mcross, 5)
for elmethod in (None, "composition", #"constant", 
                 #"coordination", "inverse_cord", 
                 "coulomb_ZiZj/d", 
                 "coulomb_1/d",):
    print elmethod
    print "knn MAE", knn_regression(mtrain, mcross, 5, elmap=elmap, elmethod=elmethod)
    print "krr MAE", krr_regression(mtrain, mcross, 50, 0.01, elmap=elmap, elmethod=elmethod)
    print "forest", sklearn_regression(mtrain, mcross, "forest", elmap=elmap, elmethod=elmethod)  
    print "svr",    sklearn_regression(mtrain, mcross, "svr", elmap=elmap, elmethod=elmethod)  


#pca_decomposition(mtrain, mcross)

# knn, forest, svr(nusvr) are the best
#for  method in (["tree", "forest", "bayesridge", "ard", "lars", "lasso", "linear", "passiveagressive", "sgd", "svr", "nusvr"]):
#    print method, sklearn_regression(mtrain, mcross, method)

# svr with rbf kernel is the best
# gamma (defaul=0) is good, degree does not matter, probablity does not matter when epsilon = 0
#for kernel in (['rbf','linear', 'poly']):
#    for epsilon in ([0, 0.01, 0.1, ]):
#        for prob in ([True, False]):
#            print kernel, epsilon, sklearn_regression(mtrain, mcross, method="svr", kernel=kernel, epsilon=epsilon, probability=prob)
#print sklearn_regression(mtrain, mcross, method="svr")









    

