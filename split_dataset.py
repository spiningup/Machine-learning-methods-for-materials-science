import numpy as np

# five fold cross validation
# mset is sorted by Eref
def get_testset(mset):
    mset = mset[:len(mset)//5*5] 
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
    
def get_partial_dataset(mset, choice=1):
    if choice == 5: return mset
    mset = mset[:len(mset)//5*5] 
    mnewset = []
    for j in range(choice):
        i = 0
        while i < len(mset):
            idx = np.random.randint(0, 5-j)
            mnewset.append(mset[i+idx])
            mset.pop(i+idx)
            i += 5 - j - 1
    return mnewset
