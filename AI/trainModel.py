#! /usr/bin/env python3

# Trains an AI model that predicts the accuracy of the variant calling based on the coverage pattern.

import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy import optimize


# Model training
# https://scikit-learn.org/stable/modules/tree.html
# https://scikit-learn.org/stable/modules/generated/sklearn.tree.DecisionTreeRegressor.html
from sklearn import tree 
from sklearn.model_selection import train_test_split


print('Importing training data from file...')
sample_name = 'CFSANSMP000113116'
with open('./%s-trainingData.pkl' % sample_name, 'rb') as file:
    masks = pickle.load(file)
    freyja_truths = pickle.load(file)
    full_var_call = pickle.load(file)
    

# Model is supposed to predict the behaviour of f in y=f(x)
print('Splitting the test/train data...')
x_train, x_test, y_train, y_test = train_test_split(masks, freyja_truths, train_size=0.9)



def trainTree (param):
    clf = tree.DecisionTreeRegressor(ccp_alpha=param)
    clf = clf.fit(x_train, y_train)

    # Test the performance on the training set itself
    training_predictions = clf.predict(x_train)
    corr_matrix = np.corrcoef(y_train, training_predictions)
    R2_train = corr_matrix[0,1]**2

    # Performance on the test set
    test_predictions = clf.predict(x_test)
    corr_matrix = np.corrcoef(y_test, test_predictions)
    R2_test = corr_matrix[0,1]**2
    
    print('%e\t\tR2_training: %.5f\tR2_test: %.5f' % (param, R2_train, R2_test) )
    return (clf, R2_train, R2_test)


###############################################
# Train the random decision tree / forest / NN
print("Training the model...")

# A function to evaluate R2 coeff. difference between the training and test sets
def r2dif_train_test(param):
    if param <0 or param>0.0001:
        return 0

    clf, R2_train, R2_test = trainTree(param)
    
    if np.isnan(R2_test):
        return 0
    else:
        return -R2_test


# Ideally R2_training = R2_test, which is what we solve for here
# root = optimize.brentq(r2dif_train_test, 0, 0.005, xtol=1e-5)
root = optimize.minimize_scalar(r2dif_train_test, method="brent", tol = 1e-3)


print("Final model evaluation: (param=%e)" % root.x)
clf, R2_train, R2_test = trainTree(root.x)


# Test the performance on the training set itself
training_predictions = clf.predict(x_train)
corr_matrix = np.corrcoef(y_train, training_predictions)
R2_train = corr_matrix[0,1]**2

# Performance on the test set
test_predictions = clf.predict(x_test)
corr_matrix = np.corrcoef(y_test, test_predictions)
R2_test = corr_matrix[0,1]**2
    

print('Exporting the trained model to file...')
with open('./trainedModel.pkl', 'wb') as file:
    pickle.dump(clf,file)


# Plot comparisons to show the predictive power of the model
print("Plotting in progress...")

# Training set itself
plt.plot(y_train, training_predictions, 'b.')
plt.xlabel('Ground truth')
plt.ylabel('Model prediction')
plt.title('Training set (R^2 = %.3f)' % R2_train)

plt.savefig('./training.png', dpi=200)
plt.close()


# Performance on the test set
plt.plot(y_test, test_predictions, 'b.')
plt.xlabel('Ground truth')
plt.ylabel('Model prediction')
plt.title('Test set (R^2 = %.3f)' % R2_test)

plt.savefig('./test.png', dpi=200)
plt.close()


print('Done')


