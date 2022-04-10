#! /usr/bin/env python3

# Trains an AI model that predicts the accuracy of the variant calling based on the coverage pattern.

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import pickle
import multiprocessing


# Model training
# https://scikit-learn.org/stable/modules/tree.html
# https://scikit-learn.org/stable/modules/generated/sklearn.tree.DecisionTreeRegressor.html
from sklearn import tree 
from sklearn.model_selection import train_test_split


def calculate_R2 (x,y):
    corr_matrix = np.corrcoef(x, y)
    R2 = corr_matrix[0,1]**2
    return R2


num_features = 29903
all_masks = np.empty((0,num_features))
all_tvs = np.empty(0)
for sample_name in ['CFSANSMP000113116']: #, 'CFSANSMP000113117', 'CFSANSMP000113118', 'CFSANSMP000113119']:
    filename = '/hpc/scratch/Tunc.Kayikcioglu/%s-trainingData.pkl' % sample_name
    print('Importing training data from %s' % filename)

    with open(filename, 'rb') as file:
        masks = pickle.load(file)
        freyja_deviations = pickle.load(file)
        full_var_call = pickle.load(file)
    
    all_masks = np.vstack((all_masks, masks)) # Coverage pattern obtained/simulated
    all_tvs = np.hstack((all_tvs, 100*np.abs(freyja_deviations-full_var_call))) # True values to predict


num_samples = len(all_tvs)


###############################################
# Take a random subset out of the available data and build a tree.
# Prune the tree to maximise the R2 on the test subset.
def makeTree (tree_id):
    print('Generating data subset for %d...' % tree_id)
    selected_entries = np.random.choice(num_samples, subset_size, replace=False)
    masks = all_masks[selected_entries]
    tvs = all_tvs[selected_entries]    

    # Model is supposed to predict the behaviour of f in y=f(x)
    # print('Splitting the test/train data...')
    x_train, x_test, y_train, y_test = train_test_split(masks, tvs, train_size=0.5)
    
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

    
    # A function to evaluate R2 coeff. difference between the training and test sets
    def r2dif_train_test(param):
        if param <0 or param>0.1:
            return 1

        clf, R2_train, R2_test = trainTree(param)
        
        if np.isnan(R2_test):
            return 1
        else:
            return R2_train-R2_test


    #Train the random decision tree / forest / NN
    # Ideally R2_test is high which is what we optimize for here
    print("Training model %d..." % tree_id)
    root = optimize.minimize_scalar(r2dif_train_test, method="brent", options={'xtol':1e-2})

    clf, R2_train, R2_test = trainTree(root.x)
    
    ### Report out the output, one thread at a time ###
    ### -------------------- ###
    critical_lock.acquire()

    print("Done with %d..." % tree_id)
    print('Alpha\t\tR2_train\t\tR2_test')
    for iter in func_evaluations:
        print('%e\t%.5f\t\t%.5f' % iter)
    print('')
    
    critical_lock.release()
    ### -------------------- ###
    
    return clf
    

# Evaluate each tree on a separate thread
critical_lock = multiprocessing.Lock()
subset_size = 2000
num_trees = 10

# Spawn mutliple processes
pool = multiprocessing.Pool()
models = pool.map(makeTree, range(num_trees))

# Wait for all processes to complete
pool.close()
pool.join()
    

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
plt.hist2d(tvs, predictions, bins=20, cmap='GnBu')
plt.xlabel('Ground truth')
plt.ylabel('Model prediction')
plt.title('Bootstrapping set (R^2 = %.5f)' % R2)

plt.savefig('./bootstrapping.png', dpi=200)
plt.close()


print('All done')


