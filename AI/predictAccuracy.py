#!/usr/bin/env python3

# Trains an AI model that predicts the accuracy of the variant calling based on the coverage pattern.

import matplotlib.pyplot as plt
import numpy as np
import sys, pickle, csv


# https://scikit-learn.org/stable/modules/tree.html
from sklearn import tree 


# Get and test input parameters
if len(sys.argv) < 3:
    raise Exception('Incorrect call to the script.')

model_filename = sys.argv[1]
sample_filenames = sys.argv[2:]
num_samples = len(sample_filenames)

# Import the trained tree from file
with open(model_filename, 'rb') as file:
    num_features = pickle.load(file)
    num_models = pickle.load(file)

    models = []
    for i in range(num_models):
        clf = pickle.load(file)
        models.append(clf)


# Import the csv file with columns pos; coverage; quality
def import_sample_features(filename):
    input_feature = np.zeros(num_features)
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter="\t")
        for row in reader:
            pos_idx = int(row[1])-1
            read_depth = row[2]
            # avg_quality = row[3]
            input_feature[pos_idx] = read_depth
    
    return input_feature


predictions = []
for filename in sample_filenames:
    feature = import_sample_features(filename)
 
    raw_predictions = []
    for clf in models:
        raw_predictions.append(float( clf.predict(feature.reshape(1,-1)) ))
    
    predictions.append(np.median(raw_predictions))


# Plot the accuracy predictions
FDAblue = (0, 124/255, 186/255)  # RGB color representation of the logo
plt.rcParams.update({'font.size': 14})

plt.bar(np.arange(1,num_samples+1), predictions, color=FDAblue)
plt.xlabel('Sample number')
plt.xticks(np.arange(1,num_samples+1, int(max(1,num_samples/10))))
plt.ylabel('Predicted deviation (percentage points)')

plt.savefig('QC_predictions.png', dpi=200)
plt.close()

