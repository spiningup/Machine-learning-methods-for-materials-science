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
    
kRR = kernel_ridge_regression(mtrain, mcross, lamda, sigma, matrixtype=1)
print kRR.run()
#kRR.choose_lamda_sigma([10, 50], [0.01, 0.001])

print knn_regression(mtrain, mcross, 5, 3, [1, 0.01, 1])




    

