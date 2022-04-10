#!/usr/bin/env python3

# Trains an AI model that predicts the accuracy of the variant calling based on the coverage pattern.

import numpy as np
import matplotlib.pyplot as plt
import pickle, os, csv, shutil
import pandas as pd


# Functions to interact with freyja
def subsample_freyja_inputs (mask, file_dir='.', file_suffix=''):
    with open('%s/freyja.depths.%s' % (file_dir,file_suffix), 'w') as infile:
        writer = csv.writer(infile, delimiter="\t")
        for i in range(num_features):
            if mask[i] == 1:
                writer.writerow(depth_file_contents[i])
            else:
                emptyRow = depth_file_contents[i]
                emptyRow[3] = 0
                writer.writerow(emptyRow)

    with open('%s/freyja.variants.tsv.%s' % (file_dir,file_suffix), 'w') as infile:
        writer = csv.writer(infile, delimiter="\t")
        writer.writerow(var_file_header)
        for (pos_idx,row) in var_file_contents:
            if mask[pos_idx] == 1:
                writer.writerow(row)
    

def parse_freyja_output (file_dir='.', file_suffix=''):
    freyja_raw = pd.read_table('%s/freyja.demix.%s' % (file_dir,file_suffix), index_col=0)
    full_var_calls = dict(eval( pd.Series(freyja_raw.loc['summarized'][0])[0].replace('inf', '0').replace('nan', '0') ))
    return full_var_calls['Omicron']


# Generate training data by bootstrapping from the imported experimental data
print('Generating features...')

if False:
    fq = "./fastq/CFSANSMP000113116_S1_L001_R1_001.fastq.gz"
    os.system("bowtie2 --no-unal --threads 4 -x ~/c-wap/covidRefSequences/wuhan -U %s -S aligned.sam" % fq)
    os.system("samtools sort aligned.sam -o sorted.bam -@ 4")
    os.system("rm aligned.sam")

    os.system("freyja variants sorted.bam --variants freyja.variants.tsv --depths freyja.depths --ref ~/c-wap/covidRefSequences/wuhan.fa")
    os.system("rm sorted.bam")
    os.system("freyja demix freyja.variants.tsv.full freyja.depths.full --output freyja.demix.full")



# Abundance estimates by Freyja for the full data set
full_var_call = parse_freyja_output(file_dir='.', file_suffix='full')


# Import the tsv and .depths file used for this computation
print('Importing original data...')
with open('freyja.depths.full', 'r') as infile:
    reader = csv.reader(infile, delimiter="\t") 
    depth_file_contents = []
    for row in reader:
        depth_file_contents.append(row)

with open('freyja.variants.tsv.full', 'r') as infile:
    reader = csv.reader(infile, delimiter="\t")
    var_file_header = next(reader)
    var_file_contents = []
    for row in reader:
        pos_idx = int(row[1])
        var_file_contents.append( (pos_idx,row) )



print('Generating subsampling mashes for the training sets...')
if os.path.exists('./freyjaOutput'):
    shutil.rmtree('./freyjaOutput')
os.mkdir('./freyjaOutput')

exit()

num_samples = 3
num_features = len(depth_file_contents) # i.e., the genome size.


masks = np.empty((num_samples,num_features))
for i in range(num_samples):
    # Random number of loci at random locations are selected to be missing
    pct_loci_missing = np.random.randint(0,101)
    loci_present = np.ones(num_features)
    loci_present[np.random.randint(0,101,num_features) < pct_loci_missing] = 0

    subsample_freyja_inputs(loci_present, file_dir='freyjaOutput', file_suffix=i)
    masks[i,:] = loci_present


print('Feeding subsampled inputs to Freyja...')
os.system("./executeFreyja.sh")


print('Parsing freyja output...')
freyja_truths = np.empty(num_samples)
for i in range(num_samples):
    freyja_truths[i] = parse_freyja_output(file_dir='freyjaOutput', file_suffix=i)

    
# https://scikit-learn.org/stable/modules/tree.html
from sklearn import tree 
from sklearn.model_selection import train_test_split

# Model is supposed to predict the behaviour of f in y=f(x)
x_train, x_test, y_train, y_test = train_test_split(masks, freyja_truths, train_size=0.9)


# Train the random decision tree / forest
print("Training the model...")
clf = tree.DecisionTreeRegressor()
clf = clf.fit(x_train, y_train)


# Export the trained tree to a file
with open('./trainedModel.pkl', 'wb') as file:
    pickle.dump(clf, file)
    pickle.dump(num_features, file)


# Plot comparisons to show the predictive power of the model
print("Self-evaluation in progress...")

# Training set itself
training_predictions = clf.predict(x_train)
plt.plot(y_train, training_predictions, 'b*')
plt.xlabel('Ground truth')
plt.ylabel('Model prediction')
plt.title('Training set')

plt.savefig('./training.png', dpi=200)
plt.close()


# Performance on the test set
test_predictions = clf.predict(x_test)
plt.plot(y_test, test_predictions, 'b*')
plt.xlabel('Ground truth')
plt.ylabel('Model prediction')
plt.title('Test set')

plt.savefig('./test.png', dpi=200)
plt.close()


print('Done')


