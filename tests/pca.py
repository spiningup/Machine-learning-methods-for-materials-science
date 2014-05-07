import os
os.chdir("../")
from read_json import read_json
from split_dataset import *
from ml import *

mset = read_json("data.json")
mtest, mset = get_testset(mset)
mtrain, mcross, mset = get_train_validation_set(mset)

Xtrain, Xcross, Etrain, Ecross = get_X(mtrain, mcross)
Xtrain, Xcross = pca_decomposition(Xtrain, Xcross, n_components=7, kernel=None)

#output
#[ 0.52311348  0.22222354  0.12777081  0.05712278  0.02736759  0.01505734
#  0.01357351] 0.986229039451
