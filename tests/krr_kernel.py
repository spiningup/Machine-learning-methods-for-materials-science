import os
from read_json import read_json, get_elements_map
from split_dataset import *
from ml import *
np.random.seed(0)

results0 = {
    "gaussian" : {
        1   : 2.230644997509021,
        2   : 0.59520440371224048,  
        5   : 0.092406396008585051, 
        10  : 0.073114579684683864, 
        50  : 0.13806670762186987, 
        100 : 0.20564880488459666, }, 
    "laplacian": {
        1   : 1.6484438973412825, 
        2   : 0.33806515558760408, 
        5   : 0.11973026611022271, 
        10  : 0.10496598445885565, 
        50  : 0.10030346100452348, 
        100 : 0.10005148232205517, }
}
    

mset = read_json("data.json")
mtest, mset = get_testset(mset)
mtrain, mcross, mset = get_train_validation_set(mset)
elmap = get_elements_map(mset)

elmethod = "composition"
for kernel in ("gaussian", 
               "laplacian",):
    for sigma in ([1,2,5,10,50,100]):
        for lamda in ([0.0001, 0.001,]):
            MAE = krr_regression(mtrain, mcross, sigma, lamda, kernel=kernel, elmap=elmap, elmethod=elmethod)
            print "krr MAE", kernel, sigma, lamda, MAE
            if lamda == 0.0001:
                assert np.abs(MAE[1] - results0[kernel][sigma]) < 1e-6
            
