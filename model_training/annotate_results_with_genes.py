#!/usr/bin/env python
"""
Given a TBA coefficients file or a TBA significance file, maps the 
motif names to gene names
"""

### imports ###
import argparse
import numpy as np
import os    
import inspect
import pandas as pd
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Given a TBA coefficients file or \
    a TBA significance file, maps the motif names to gene names')
    parser.add_argument("result_path",
        help="path to a TBA coefficients or significance file",
        type = str)
    parser.add_argument("output_path",
        help="file path where output should be written",
        type=str)
    
    # parse arguments
    args = parser.parse_args()

    result_path = args.result_path
    output_path = args.output_path

    metadata_path= os.path.dirname(
        os.path.abspath(
        inspect.getfile(
        inspect.currentframe()))).replace('model_training', 'motif_metadata.tsv')
    
    metadata_frame = pd.read_csv(metadata_path, sep='\t')
    motifName_gene_dict = dict(zip(metadata_frame['Name'].values, 
        metadata_frame['Gene'].values))
   
    results_frame = pd.read_csv(result_path, sep='\t', index_col=0)
    columns = list(results_frame.columns.values)
    results_frame['genes'] = [motifName_gene_dict[x] for x in results_frame.index.values]
    columns.insert(0, 'genes')
    results_frame = results_frame[columns]
    results_frame.to_csv(output_path, sep='\t')
