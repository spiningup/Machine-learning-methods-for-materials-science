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
print "krr MAE", krr_regression(mtrain, mcross, 50, 0.01, elmap=elmap, elmethod=elmethod)    

