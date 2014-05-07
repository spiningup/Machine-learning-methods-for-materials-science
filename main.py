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

elmethod = "composition"
#for kernel in ("gaussian", 
#               "laplacian",):
#    for sigma in ([1,2,5,10,50,100]):
#        for lamda in ([0.0001, 0.001, 0.01, ]):
#            print "krr MAE", kernel, sigma, lamda, krr_regression(mtrain, mcross, sigma, lamda, kernel=kernel, elmap=elmap, elmethod=elmethod)    


sigma = 10 ; lamda = 0.0001 ; kernel = "gaussian"
krr_regression(mtrain, mcross, sigma, lamda, kernel=kernel, elmap=elmap, elmethod=elmethod)    
