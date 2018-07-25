#!/usr/bin/env python
"""
Given a fasta file, and a set of PWMs in Homer format calculates the top motif
match for each motif for each sequence:
"""

### imports ###
import argparse
import os
import numpy as np
import time
import multiprocessing
import pickle
import pandas as pd
import Bio
from Bio import motifs
from Bio import SeqIO

### functions ###
def read_jaspar_motif_file(motifPath, pseudocount):
    ''' 
    reads jaspar motif file
    inputs: path to a jaspar motif file
    outputs: a tuple representing a motif
    '''
    with open(motifPath) as f:
        m = motifs.read(f, 'jaspar')
        default_pseudocount = motifs.jaspar.calculate_pseudocounts(m)
        scaled_pseudocount = pseudocount/0.01 * default_pseudocount['A']
        m.pseudocounts = int(scaled_pseudocount)
    return (m.name, m)

def read_fasta(file_path):
    '''
    reads in a fasta file and returns a list of sequence ids and a list of sequences
    inputs: filepath - path to a fasta file
    outputs: sequence_list - a list of sequences
             id_list - a list of ids
    '''
    # read in sequences
    id_list = []
    sequence_list = []

    alphabet = Bio.Seq.IUPAC.Alphabet.IUPAC.IUPACUnambiguousDNA()
    for seq_record in SeqIO.parse(file_path, "fasta"):
        seq_record.seq.alphabet = alphabet

        id_list.append(seq_record.id)
        sequence_list.append(seq_record.seq)
    return sequence_list, id_list

def calculate_all_motif_scores_async(sequence_list,
                                 pssm,
                                 motif_name,
                                 output_path
                                 ):
    '''
    identifies the highest scoring match to pwm for each sequence
    inputs: pwm - a numpy array representing a pwm
            sequence_array_list - a list of numpy array representations of a sequence with:
                                 A = [1, 0, 0, 0]
                                 C = [0, 1, 0, 0]
                                 G = [0, 0, 1, 0]
                                 T = [0, 0, 0, 1]
    outputs: top_scores - a list of the best motif scores in each sequence
    '''

    fwd_pssm = pssm
    rev_pssm = fwd_pssm.reverse_complement()

    all_fwd_scores = [] 
    all_rev_scores = []

    start = time.time()

    for seq in sequence_list:
        fwd_scores = fwd_pssm.calculate(seq)
        rev_scores = rev_pssm.calculate(seq)

        all_fwd_scores.append(fwd_scores)
        all_rev_scores.append(rev_scores)

    pickle.dump(all_fwd_scores,
                open(output_path+'/'+motif_name+'_fwd.pickle','wb'))
    pickle.dump(all_rev_scores,
                open(output_path+'/'+motif_name+'_rev.pickle','wb'))

    end = time.time()
    print(motif_name, 'calculation time:', end-start)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='given a fasta file and a list \
                                     and a list of motif files, calculates the \
                                     the best match to each motif for each \
                                     sequence' )
    parser.add_argument("fasta_path",
        help="path to a fasta_file containing sequences to score",
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
        default=4)
    parser.add_argument("-pseudocount", 
        help="pseudocount for calculating motif scores",
        type=float,
        default=0.01)

    # parse arguments
    args = parser.parse_args()

    fasta_path = args.fasta_path
    output_path = args.output_path
    motif_files = args.motif_files
    num_processors = args.num_procs
    pseudocount = args.pseudocount
    
    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    # read in motif files
    all_motifs = []
    for m in motif_files:
        motif = read_jaspar_motif_file(m, pseudocount)
        all_motifs.append(motif)
    # sort motifs by name
    all_motifs.sort(key=lambda x:x[0])

    start = time.time()
    
    sequence_list, id_list = read_fasta(fasta_path)

    # convert strings to arrays
    #sequence_array_list = convert_sequences_to_array(sequence_list)
    
    # calculate motif scores
    pool = multiprocessing.Pool(processes=num_processors)
    manager = multiprocessing.Manager()
       
    for motif in all_motifs:
        motif_name = motif[0]
        pssm = motif[1].pssm
        
        pool.apply_async(
        calculate_all_motif_scores_async ,args=(sequence_list, 
                                                pssm, 
                                                motif_name, 
                                                output_path
                                           )
                        )
    pool.close()
    pool.join()
    
    end = time.time()
    print('total time', end-start)
