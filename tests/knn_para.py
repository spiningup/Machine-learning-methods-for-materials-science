import os
from read_json import read_json
from split_dataset import *
from ml import *
np.random.seed(0)

mset = read_json("data.json")
mtest, mset = get_testset(mset)
mtrain, mcross, mset = get_train_validation_set(mset)

# test knn, kernel=None, scaling=1, weights=distance, metric="minkowsk" are the best options
for kernel in ([None, 'rbf', 'poly','cosine']):
    for scaling in ([1, 2, 3]):
        for weights in (["distance", "uniform"]):
            crosserror = knn_regression(mtrain, mcross, 5, kernel=kernel, scaling=scaling, weights=weights)
            print kernel, scaling, weights, crosserror
            
            if kernel == "rbf" and scaling == 1 and weights == "distance": assert np.abs(crosserror - 0.276098175986) < 1e-6

# distance metric
for metric in (["euclidean", "manhattan", "chebyshev", "minkowski",]):
    print metric, knn_regression(mtrain, mcross, 5, metric=metric)

# select features
print "knn MAE", knn_regression(mtrain, mcross, 5, selectf=True) 


# full output 
#Size of dataset :  1155
#None 1 distance 0.276098175986
#None 1 uniform 0.291063464644
#None 2 distance 0.300953090769
#None 2 uniform 0.31007536633
#None 3 distance 0.528338033874
#None 3 uniform 0.529301697634
#rbf 1 distance 0.276098175986
#rbf 1 uniform 0.291063464644
#rbf 2 distance 0.300953090769
#rbf 2 uniform 0.31007536633
#rbf 3 distance 0.528338033874
#rbf 3 uniform 0.529301697634
#poly 1 distance 0.276098175986
#poly 1 uniform 0.291063464644
#poly 2 distance 0.300953090769
#poly 2 uniform 0.31007536633
#poly 3 distance 0.528338033874
#poly 3 uniform 0.529301697634
#cosine 1 distance 0.276098175986
#cosine 1 uniform 0.291063464644
#cosine 2 distance 0.300953090769
#cosine 2 uniform 0.31007536633
#cosine 3 distance 0.528338033874
#cosine 3 uniform 0.529301697634
#euclidean 0.276098175986
#manhattan 0.275781468692
#chebyshev 0.28520854345
#minkowski 0.276098175986
#knn MAE found
#[0, 1, 2, 3, 4, 5, 6, 7, 8] 0.273827896939
#found
#[0, 1, 2, 3, 4, 6, 7, 8] 0.273113703553
#not found
#[0, 1, 2, 3, 4, 6, 7, 8] 0.273113703553
#0.273113703553
