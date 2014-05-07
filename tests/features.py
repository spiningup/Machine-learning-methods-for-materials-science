import os
os.chdir("../")
from read_json import read_json, get_elements_map
from split_dataset import *
from ml import *

results0 = {
    None: {
        "knn": 0.276098175986, 
        "krr": 0.24158817150801776, 
        "forest" : 0.277859153956, 
        "svr" : 0.277415910627, }, 
    "composition" : {
        "knn": 0.202293656886,
        "krr": 0.11052597542220206, 
        "forest": 0.186148071251,
        "svr": 0.135622957115, }, 
    "coulomb_ZiZj/d" : {
        "knn": 0.296029527716, 
        "krr":  0.20612285173126219, 
        "forest": 0.206198925628, 
        "svr": 0.245838741856, }, 
    "coulomb_1/d" : {
        "knn": 0.320368738508, 
        "krr": 0.21328211175452691, 
        "forest": 0.20198798341, 
        "svr": 0.259278560667, }, 
}


mset = read_json("data.json")
mtest, mset = get_testset(mset)
mtrain, mcross, mset = get_train_validation_set(mset)
elmap = get_elements_map(mset)

results = {}
for elmethod in (None, "composition", #"constant", 
                 #"coordination", "inverse_cord", 
                 "coulomb_ZiZj/d", 
                 "coulomb_1/d",):
    results[elmethod] = {}
    results[elmethod]["knn"] = knn_regression(mtrain, mcross, 5, elmap=elmap, elmethod=elmethod)
    results[elmethod]["krr"] = krr_regression(mtrain, mcross, 50, 0.01, elmap=elmap, elmethod=elmethod)[1]
    results[elmethod]["forest"] = sklearn_regression(mtrain, mcross, "forest", elmap=elmap, elmethod=elmethod)  
    results[elmethod]["svr"] = sklearn_regression(mtrain, mcross, "svr", elmap=elmap, elmethod=elmethod)  

    for key in results[elmethod].keys():
        assert np.abs(results[elmethod][key] - results0[elmethod][key]) < 1e-6

