import os
from read_json import read_json
from split_dataset import *
from ml import *
np.random.seed(0)

mset = read_json("data.json")
mtest, mset = get_testset(mset)
mtrain, mcross, mset = get_train_validation_set(mset)

# svr with rbf kernel is the best
# gamma (defaul=0) is good, degree does not matter, probablity does not matter when epsilon = 0
for kernel in (['rbf','linear', 'poly']):
    for epsilon in ([0, 0.01, 0.1, ]):
        for prob in ([True, False]):
            print kernel, epsilon, sklearn_regression(mtrain, mcross, method="svr", kernel=kernel, epsilon=epsilon, probability=prob)

assert np.abs(sklearn_regression(mtrain, mcross, method="svr") - 0.277415910627) < 1e-6

# full output
#rbf 0 0.277924113858
#rbf 0 0.277924113858
#rbf 0.01 0.277677090125
#rbf 0.01 0.277677090125
#rbf 0.1 0.277415910627
#rbf 0.1 0.277415910627
#linear 0 0.374431723306
#linear 0 0.374431723306
#linear 0.01 0.373983540223
#linear 0.01 0.373983540223
#linear 0.1 0.373714471369
#linear 0.1 0.373714471369
#poly 0 0.436654104565
#poly 0 0.436654104565
#poly 0.01 0.438642394784
#poly 0.01 0.438642394784
#poly 0.1 0.438529691286
#poly 0.1 0.438529691286
#0.277415910627
