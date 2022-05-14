#! /usr/bin/env python3

# Trains an AI model that predicts the accuracy of the variant calling based on the coverage pattern.

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import pickle
from multiprocessing import Pool, Lock, shared_memory


# Model training
# https://scikit-learn.org/stable/modules/tree.html
# https://scikit-learn.org/stable/modules/generated/sklearn.tree.DecisionTreeRegressor.html
from sklearn import tree 
#from sklearn.model_selection import train_test_split


num_features = 29903
all_masks = np.empty((0,num_features))
all_tvs = np.empty(0)
for sample_name in ['CFSANSMP000113116', 'CFSANSMP000113117', 'CFSANSMP000113118', 'CFSANSMP000113119']:
    filename = '/hpc/scratch/Tunc.Kayikcioglu/%s-trainingData.pkl' % sample_name
    print('Importing training data from %s' % filename)

    with open(filename, 'rb') as file:
        masks = pickle.load(file)
        freyja_deviations = pickle.load(file)
        full_var_call = pickle.load(file)
    
    all_masks = np.vstack((all_masks, masks)) # Coverage pattern obtained/simulated
    all_tvs = np.hstack((all_tvs, 100*np.abs(freyja_deviations-full_var_call))) # True values to predict


# Generate a shared memory to store the input DB so that the data does not get copied
# once per each thread
def make_shared_memory (x):
    shm = shared_memory.SharedMemory(create=True, size=x.nbytes)
    x_shm = np.ndarray(x.shape, dtype=x.dtype, buffer=shm.buf)
    x_shm[:] = x[:]
    return (x_shm,shm)


print('Setting shared memory...')
all_masks, shm_masks = make_shared_memory(all_masks)
all_tvs, shm_tvs = make_shared_memory(all_tvs)
num_samples = len(all_tvs)


def calculate_R2 (x,y):
    corr_matrix = np.corrcoef(x, y)
    R2 = corr_matrix[0,1]**2
    return R2


###############################################
# Take a random subset out of the available data and build a tree.
# Prune the tree to maximise the R2 on the test subset.
def makeTree (tree_id):
    print('Generating data subset for %d...' % tree_id)
    np.random.seed(tree_id)
    train_selected_entries = np.random.choice(num_samples, train_size, replace=False)
    x_train = all_masks[train_selected_entries]
    y_train = all_tvs[train_selected_entries]

    test_selected_entries = np.random.choice(num_samples, test_size, replace=False)
    x_test = all_masks[test_selected_entries]
    y_test = all_tvs[test_selected_entries]

    # Split the data into test/training sets
    # x_train, x_test, y_train, y_test = train_test_split(masks, tvs, train_size=0.8) 
    

    # Model is supposed to predict the behaviour of f in y=f(x)    
    func_evaluations = []
    def trainTree (param):
        clf = tree.DecisionTreeRegressor(ccp_alpha=param, random_state=0)
        clf = clf.fit(x_train, y_train)

        # Test the performance on the training set itself
        training_predictions = clf.predict(x_train)
        R2_train = calculate_R2(y_train, training_predictions)

        # Performance on the test set
        test_predictions = clf.predict(x_test)
        R2_test = calculate_R2(y_test, test_predictions)

        func_evaluations.append( (param, R2_train, R2_test) )
        return (clf, R2_train, R2_test)

    
    # A function to evaluate R2 coeff. on the test data set to avoid over-learning.
    def r2dif_train_test(param):
        if param <0 or param>0.2:
            return 1

        clf, R2_train, R2_test = trainTree(param)
        
        if np.isnan(R2_test):
            return 1
        else:
            return -R2_test


    #Train the random decision tree / forest / NN
    # Ideally R2_test is high which is what we optimize for here
    print("Training model %d..." % tree_id)
    root = optimize.minimize_scalar(r2dif_train_test, method="brent", options={'xtol':1e-2})
    final_param = root.x

    clf, R2_train, R2_test = trainTree(final_param)
    func_evaluations.append((final_param, R2_train, R2_test))
    
    ### Report out the output, one thread at a time ###
    ### -------------------- ###
    critical_lock.acquire()

    print('')
    print("Done with %d..." % tree_id)
    print('Alpha\t\tR2_train\tR2_test')
    for iter in func_evaluations:
        print('%e\t%.5f\t\t%.5f' % iter)
    print('')
    
    critical_lock.release()
    ### -------------------- ###
    
    return (clf, R2_test)
    

# Evaluate each tree on a separate thread
critical_lock = Lock()
train_size = 50000
test_size = 10000
num_trees = 100

# Spawn multiple processes
pool = Pool(processes=20, maxtasksperchild=1)
models_R2 = pool.map(makeTree, range(num_trees))

# Wait for all processes to complete
pool.close()
pool.join()

# Sort the models w.r.t. R2_test and keep the best (i.e. highest R2) ones only
models_R2.sort(key = lambda x: x[1], reverse=True)
minimum_R2 = 0.8*models_R2[0][1]
models = [ x[0] for x in models_R2 if x[1] >= minimum_R2 ]
print('Of %d models, %d were retained.' % (len(models_R2), len(models)) )
print('    R2_test max: %.5f' % models_R2[0][1])
print('    R2_test min: %.5f' % models_R2[-1][1])
print('    R2_test limit: %.5f' % minimum_R2)
    

print('Exporting the trained models to file...')
with open('./trainedModel.pkl', 'wb') as file:
    pickle.dump(num_features,file)
    pickle.dump(len(models),file)
    for clf in models:
        pickle.dump(clf,file)
        

def multi_model_prediction(data):
    predictions = []
    for clf in models:
        predictions.append(clf.predict(data))
    
    return np.median(predictions)


print('Assessing the final quality by bootstrapping...')
test_size = 100000
selected_entries = np.random.choice(num_samples, test_size, replace=False)
predictions = [ multi_model_prediction(all_masks[x,:].reshape(1,-1)) for x in selected_entries ]
tvs = all_tvs[selected_entries]
R2 = calculate_R2(predictions,tvs)
print('R2 across the final test case: %.5f' % R2)

# Plot comparisons to show the predictive power of the model
print("Plotting in progress...")
plt.hist2d(tvs, predictions, bins=np.arange(0,101,5), cmap='GnBu')
plt.xlabel('Ground truth')
plt.xlim(0,100)
plt.ylabel('Model prediction')
plt.ylim(0,100)
plt.title('Bootstrapping set (R^2 = %.5f)' % R2)

plt.savefig('./bootstrapping.png', dpi=200)
plt.close()


print('All done')
shm_masks.close()
shm_masks.unlink()

shm_tvs.close()
shm_tvs.unlink()

