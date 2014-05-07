import os
from read_json import read_json
from split_dataset import *
from ml import *
np.random.seed(0)

mset = read_json("data.json")
mtest, mset = get_testset(mset)
mtrain, mcross, mset = get_train_validation_set(mset)

# knn, forest, svr(nusvr) are the best
result = {}
for  method in (["tree", "forest", "bayesridge", "ard", "lars", "lasso", "linear", "passiveagressive", "sgd", "svr", "nusvr"]):
    result[method] =  sklearn_regression(mtrain, mcross, method)

# target result
result0 = {
    "tree"       :  0.322621783437,     
    "forest"     :  0.27409744122 ,
    "bayesridge" :  0.374720207598,
    "ard"        :  0.375426020771,
    "lars"       :  1.02694916818 ,
    "lasso"      :  0.924103892157,
    "linear"     :  0.375011797537,
    "passiveagressive" :  1.07601888579, 
    "sgd"        :  0.445285999347,
    "svr"        :  0.277415910627,
    "nusvr"      :  0.282320263535, 
}


for key, value in result.items():
    assert np.abs(value - result0[key]) < 1e-6
