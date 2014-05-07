import numpy as np
from read_json import read_json, get_elements_map
from split_dataset import *
from ml import *

mset = read_json("data.json")
#mtest, mset = get_testset(mset)
#mtrain, mcross, mset = get_train_validation_set(mset)
mcross, mtrain = get_testset(mset)
elmap = get_elements_map(mset)

elmethod = "composition"
sigma = 10 ; lamda = 0.0001 ; kernel = "gaussian"
print krr_regression(mtrain, mcross, sigma, lamda, kernel=kernel, elmap=elmap, elmethod=elmethod)    
