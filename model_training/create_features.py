#!/usr/bin/env python
"""
Given a set of negative  sequences and positive sequences (in FASTA format) 
as well as a list of motifs, calculates motif scores 
and sequence labels suitable for training classifier
"""

### imports ###
import argparse
import numpy as np
import os
import time
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Given a set of negative \
                sequences and positive sequences (in FASTA format) as well \
                as a list of motifs, calculates motif scores \
                and sequence labels suitable for training classifier' )
    parser.add_argument("positive_sequences_path",
        help="path to a fasta_file containing positive sequences to score",
        type = str)
    parser.add_argument("negative_sequences_path",
        help="path to a fasta_file containing negative sequences to score",
        type = str)
    parser.add_argument("output_path",
        help="directory where output file should be written",
        default="./", type=str)
    parser.add_argument("motif_files",
        help="list of motif files",
        type=str,
        nargs="+")
    parser.add_argument("-num_procs", 
        help="number of processor cores to use",
        type=int,
        default=8)
    parser.add_argument("-pseudocount", 
        help="pseudocount for calculating motif scores",
        type=float,
        default=0.01)

    # parse arguments
    args = parser.parse_args()

    positive_sequences_path = args.positive_sequences_path
    negative_sequences_path = args.negative_sequences_path
    output_path = args.output_path
    motif_files = args.motif_files
    num_processors = args.num_procs
    pseudocount = args.pseudocount

    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    start = time.time()
    script_path = os.path.dirname(__file__)
    
    print('computing features for positive sequences')
    os.system(' '.join(['python', 
                        script_path+'/calculate_motif_scores_biopython.py', 
                        positive_sequences_path, 
                        output_path, ' '.join(motif_files).replace(')', '\)').replace('(', '\('),
                        '-num_procs', str(num_processors),
                        '-pseudocount', str(pseudocount), 
                        ]))

    print('computing features for negative sequences')
    os.system(' '.join(['python', 
                        script_path+'/calculate_motif_scores_biopython.py', 
                        negative_sequences_path, 
                        output_path, ' '.join(motif_files).replace(')', '\)').replace('(', '\('),
                        '-num_procs', str(num_processors),
                        '-pseudocount', str(pseudocount), 
                        ]))

    positive_name_root = positive_sequences_path.split('/')[-1].split('.')[0]
    positive_score_path = output_path + '/' + positive_name_root + '_motif_scores.tsv'

    negative_name_root = negative_sequences_path.split('/')[-1].split('.')[0]
    negative_score_path = output_path + '/' + negative_name_root + '_motif_scores.tsv'

    positive_score_frame = pd.read_csv(positive_score_path,sep='\t',index_col=0)
    negative_score_frame = pd.read_csv(negative_score_path,sep='\t',index_col=0)

    # concatenate data frames
    combined_frame = pd.concat([positive_score_frame, negative_score_frame])

    feature_out_path = output_path +  '/combined_features.tsv'
    combined_frame.to_csv(feature_out_path, sep='\t')
    
    # create labels
    print('creating labels')
    label_path = output_path + '/labels.txt'
    label_file = open(label_path, 'w')
    for ind in positive_score_frame.index.values:
        label_file.write(ind + '\t1\n')
    for ind in negative_score_frame.index.values:
        label_file.write(ind + '\t0\n')
    label_file.close()

    end = time.time()
