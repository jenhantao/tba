#!/usr/bin/env python
"""
Given a fasta file, and a set of PWMs in Homer format calculates the top motif
match for each motif for each sequence:
"""

### imports ###
import argparse
import numpy as np
import time
import multiprocessing
import pandas as pd
import os
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


def calculate_top_motif_matches_async(sequence_list, 
                                      pssm, 
                                      motif_name, 
                                      motif_score_dict, 
                                      motif_start_dict,
                                     ):
    '''
    identifies the highest scoring match to pssm for each sequence
    inputs:    
    outputs: top_scores - a list of the best motif scores in each sequence
             top_starts - a list of the start position of the best motif match in each sequence
    '''
    start = time.time()
    
    fwd_pssm = pssm
    rev_pssm = fwd_pssm.reverse_complement()
    
    top_scores = [] # motif score of best match for each sequence
    top_starts = [] # start position of best match for each sequence
    
    # calculate scores for each motif at each position
    for seq in sequence_list:
        fwd_scores = fwd_pssm.calculate(seq)
        rev_scores = rev_pssm.calculate(seq)
        
        max_fwd_score = np.max(fwd_scores)
        max_rev_score = np.max(rev_scores)
        
        if max_fwd_score > 0 or max_rev_score > 0:
            if max_fwd_score > max_rev_score:
                top_scores.append(max_fwd_score)
                top_starts.append(str(np.argmax(fwd_scores)) + ' +')
            else:
                top_scores.append(max_rev_score)
                top_starts.append(str(np.argmax(rev_scores)) + ' -') 
        else:
            top_scores.append(0)
            top_starts.append('-1 ?')

    motif_score_dict[motif_name] = top_scores
    motif_start_dict[motif_name] = top_starts
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
    motif_score_dict= manager.dict() # {motif_name:[scores]}
    motif_start_dict = manager.dict() # {motif_name:[motif_start_pos]}
       
    for motif in all_motifs:
        pssm = motif[1].pssm
        motif_name = motif[0]
        
        pool.apply_async(
        calculate_top_motif_matches_async,args=(sequence_list, 
                                                pssm, 
                                                motif_name, 
                                                motif_score_dict, 
                                                motif_start_dict
                                           )
                        )
    pool.close()
    pool.join()

    motif_score_dict = dict(motif_score_dict)
    motif_start_dict = dict(motif_start_dict)
    motif_score_frame = pd.DataFrame(motif_score_dict, 
                                     index = id_list)
    motif_start_frame = pd.DataFrame(motif_start_dict, 
                                     index = id_list)
    name_root = fasta_path.split('/')[-1].split('.')[0]
    score_tsv_path = output_path + '/' + name_root + '_motif_scores.tsv'
    start_tsv_path = output_path + '/' + name_root + '_motif_starts.tsv'
    motif_score_frame.to_csv(score_tsv_path, sep='\t')
    motif_start_frame.to_csv(start_tsv_path, sep='\t')
    
    end = time.time()
    print('total time', end-start)
