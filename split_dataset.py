import numpy as np

# five fold cross validation
# mset is sorted by Eref

np.random.seed(0)

def get_testset(mset):
    mtest = []
    i = 0
    while i < len(mset):
        idx = np.random.randint(0, 5)
        mtest.append(mset[i+idx])
        mset.pop(i+idx)
        i += 4
        
    return mtest, mset

def get_train_validation_set(mset):
    mcross = []
    mtrain = mset[:]
    i = 0
    while i < len(mtrain):
        idx = np.random.randint(0, 4)
        mcross.append(mtrain[i+idx])
        mtrain.pop(i+idx)
        i += 3
        
    return mtrain, mcross, mset # mset is a sum of mtrain and mcross
    
