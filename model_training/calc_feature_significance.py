#!/usr/bin/env python
"""
Given standardized motif features and matching labels trains a classifier and 
returns performance metrics and model coefficients
"""

### imports ###
import argparse
import numpy as np
import os
import time
import pandas as pd
from sklearn import preprocessing
import sklearn
from sklearn import linear_model
from sklearn.model_selection import train_test_split
import scipy
from joblib import Parallel, delayed

### functions ###

def calc_model_log_likelihood(probas, labels):
    log_likelihood = 0
    Y = labels.astype(float)
    for i in range(len(Y)):
        p_true = probas[i][1]
        p_false = probas[i][0]
        y = Y[i]
        prod = ((p_true)**y) * ((p_false) ** (1-y))
        log_prod = np.log(prod)
        log_likelihood += log_prod
    return log_likelihood

def calc_feature_pvals(
    features,
    labels,
    test_size=0.2,
    num_iterations=5,
    num_procs=4
    ):
    pvals = []
    num_motifs = features.shape[1]
    # split data into training and test sets
    scaler = preprocessing.StandardScaler()

    # standardize features
    standardized_features = pd.DataFrame(scaler.fit_transform(features))
    standardized_features.columns = features.columns.values
    standardized_features.index = features.index.values

    for i in range(num_iterations):
        training_features, test_features, training_labels, test_labels = train_test_split(
            features, labels, test_size = test_size)
        
        # standardize training features
        standardized_training_features = pd.DataFrame(scaler.fit_transform(training_features))
        standardized_training_features.columns = training_features.columns.values
        standardized_training_features.index = training_features.index.values
            
        #  Train affinity classifier
        classifier = sklearn.linear_model.LogisticRegression(penalty='l1', 
            solver='liblinear') 

        classifier.fit(standardized_training_features, training_labels)
        # score predictions

        probas = classifier.predict_proba(standardized_features) # [[p_false, p_true]...] 
        overall_log_likelihood = calc_model_log_likelihood(probas, labels)
        iter_pvals = []

        iter_pvals = Parallel(n_jobs=num_procs)(
            delayed(train_perturbed_classifier)(
                standardized_features, 
                labels, 
                standardized_training_features, 
                training_labels, 
                motif_to_drop,
                overall_log_likelihood) for motif_to_drop in features.columns.values)

        pvals.append(iter_pvals)

    return pvals

           

def train_perturbed_classifier(features, 
    labels, 
    training_features, 
    training_labels, 
    motif_to_drop, 
    overall_log_likelihood):

    start = time.time()
    current_features = features.drop(motif_to_drop, axis=1, inplace=False)
    current_training_features = training_features.drop(motif_to_drop, axis=1, inplace=False)
    current_classifier = sklearn.linear_model.LogisticRegression(penalty='l1', 
        solver='liblinear')
    current_classifier.fit(current_training_features, training_labels)

    current_probas = current_classifier.predict_proba(current_features)
    current_log_likelihood = calc_model_log_likelihood(current_probas, labels)
    stat = -2*(current_log_likelihood - overall_log_likelihood)
    p = scipy.stats.chisqprob(stat, df=1)
    print('tested', motif_to_drop, time.time() - start)
    return p

    
def read_labels(label_path):
    '''
    reads label files created by create_features.py and returns a pandas Series representation
    '''
    indices = []
    vals = []
    with open(label_path) as f:
        data = f.readlines()
    for line in data:
        tokens = line.strip().split()
        indices.append(tokens[0])
        if tokens[1] == '1':
            vals.append(True)
        else:
            vals.append(False)
    to_return = pd.Series(vals, index=indices)
    return to_return

def write_test_results(features, pvals, output_path):
    pval_dict = dict(zip(range(len(pvals)), pvals))
    pval_frame = pd.DataFrame(data = pval_dict, index = features.columns.values)
    pval_frame.to_csv(output_path + '/significance.tsv', sep = '\t')
    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Performs an in silico mutagenesis test to assign significance to each motif' )
    parser.add_argument("feature_path",
        help="path to a standardized feature set created by create_features.py",
        type = str)
    parser.add_argument("label_path",
        help="path to a fasta_file containing negative sequences to score",
        type = str)
    parser.add_argument("output_path",
        help="directory where output file should be written",
        default="./", type=str)
    parser.add_argument("-num_iterations", 
        help="number of iterations to train classifier",
        type=int,
        default=5)
    parser.add_argument("-test_fraction", 
        help="fraction of data to use for testing classifier",
        type=float,
        default=0.2)
    parser.add_argument("-num_procs", 
        help="number of processors to use",
        type=int,
        default=4)

    # parse arguments
    args = parser.parse_args()

    feature_path = args.feature_path
    label_path = args.label_path
    output_path = args.output_path
    num_iterations = args.num_iterations
    test_fraction = args.test_fraction
    num_procs = args.num_procs

    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    # read in features
    print('reading features and labels')
    feature_frame = pd.read_csv(feature_path, sep='\t', index_col=0)

    # read in labels
    labels = read_labels(label_path)

    print('testing feature significance for', feature_path)
    pvals = calc_feature_pvals(feature_frame,
                            labels,
                            num_iterations=num_iterations,
                            num_procs = num_procs,
                            test_size=test_fraction)
    print('writing results')
    write_test_results(feature_frame, pvals, output_path) 
