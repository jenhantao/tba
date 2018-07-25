#!/usr/bin/env python
"""
Given a set of genomic coordinates in BED format:
chr start end
...

calculates a GC matched set of genomic coordinates.
Only first 3 columns of input file will be used- all other columns are ignored.
"""

### imports ###
import os
import sys
import numpy as np
import argparse
import inspect


def read_target_positions(file_path, filter_chromosomes):
    """
    reads a bed file and returns a list of tuples containing genomic coordinates
    """
    if filter_chromosomes==None:
        filter_chromosomes = []
    else:
        print('filtering out: ' + ' '.join(filter_chromosomes))
    with open(file_path) as f:
        data = f.readlines()
    filter_chromosomes = set(filter_chromosomes)
    positions = []
    if data[0].strip()[0] == '#':
        data = data[1:]
    for line in data:
        tokens = line.strip().split()
        chrom = tokens[0]
        start = int(tokens[1])
        end = int(tokens[2])
        #name = tokens[3]
        if not chrom in filter_chromosomes:
            if not 'chrUn' in chrom and not 'random' in chrom and not 'alt' in chrom:
                positions.append([chrom, start, end])
    return positions 

def calc_gc_content(sequence):
    '''
    sequence - a string, representing a DNA sequence in upper case
    returns the GC content of sequence
    '''
    C_count = sequence.count('C')
    G_count = sequence.count('G')
    GC_count = C_count + G_count
    GC_content = GC_count/len(sequence)
    return GC_content

def get_random_background(target_positions,
                          size_ratio,
                          num_bins = 10,
                          n_threshold = 0.5,
                          genome = 'mm10',
                          filter_chromosomes = ['chrM', 'chrY']
                          ):
    """
    target_sequences: 2D numpy array, list of genomic coordinates for target 
                      sequences [[chr,start,end],...]
    size_ratio: float, ratio of background sequences to target sequences
    num_bins: int, number of GC bins
    n_threshold: proportion of background sequences that can be N
    genome: genome from which to draw background sequences
    """
    
    ###load genome into memory
    
    # index target positions
    # {chr:[]}, value is chromosome length boolean array
    # largest chromosome has 200 million bps 
    script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    genome_path = script_path + '/' + genome + '/'

    chromosomes = [x.split('.')[0] for x in os.listdir(genome_path)]
    chromosomes = [chrom for chrom in chromosomes if not 'chrUn' in chrom and not 'random' in chrom and not 'alt' in chrom]
    filter_chromosomes = set(filter_chromosomes)
    chromosomes = [chrom for chrom in chromosomes if not chrom in filter_chromosomes] 

    chrom_size_dict = {}
    chrom_seq_dict = {}

    print('reading genome', genome)
    for chrom in chromosomes:
        with open(genome_path + chrom + '.fa') as f:
            data = f.readlines()
        seq = ''.join(x.upper().strip() for x in data[1:])
        size = len(seq)
        chrom_size_dict[chrom] = size
        chrom_seq_dict[chrom] = seq

    print('done reading genome')

    ### initialize target_chr_position_dict using target positions
    target_chr_position_dict = {x:np.zeros(200000000) for x in chromosomes} 
    ### retrieve target sequences and calculate GC content and mean length
    target_length_count = 0
    filtered_target_positions = []
    for pos in target_positions:
        chrom = pos[0]        
        start = int(pos[1])
        end = int(pos[2])
        # use 0 indexing of position, versus 1 indexing used in fasta
        if chrom in chrom_seq_dict:
            seq = chrom_seq_dict[chrom][start:(end)]
            target_chr_position_dict[chrom][start-1:end] = 1         
            if len(seq) > 0:
                gc_content = calc_gc_content(seq)
                pos.append(seq)
                pos.append(gc_content)
                filtered_target_positions.append(pos)
                target_length_count += len(seq)
            else:
                print(chrom, start, end, 'not found')
        else:
            print(chrom, start, end, 'not found')
            

    # average length of target sequences
    mean_target_length = target_length_count/len(filtered_target_positions)     
    mean_target_length = int(np.floor(mean_target_length))

    # sort target_positions by gc_content and bin according to GC content
    sorted_target_positions = sorted(filtered_target_positions, key=lambda x:x[-1])
    sorted_target_positions = np.array(sorted_target_positions)

    target_position_bins = np.array_split(sorted_target_positions, num_bins)

    min_gc = float(sorted_target_positions[0][-1])
    max_gc = float(sorted_target_positions[-1][-1])

    gc_threshold = (max_gc - min_gc)/(num_bins*2)

    background_positions = [] 

    for target_pos_bin in target_position_bins:
        current_random_pos = get_random_positions_with_gc(target_pos_bin,
                                                         size_ratio,
                                                         gc_threshold,
                                                         n_threshold,
                                                         chrom_seq_dict,
                                                         chrom_size_dict,
                                                         target_chr_position_dict,
                                                         mean_target_length)
        background_positions = background_positions + current_random_pos

    return background_positions

def get_random_positions_with_gc(target_positions, 
                                 size_ratio, 
                                 tolerance, 
                                 n_threshold,
                                 chrom_seq_dict,
                                 chrom_size_dict,
                                 target_chr_position_dict,
                                 mean_target_length
                                 ):
    """
    target_positions: 2D numpy array, list of genomic coordinates for target 
                      sequences [[chr,start,end, seq, gc_content],...]
    size_ratio: float, ratio of background sequences to target sequences
    tolerance: float, max difference in GC content between target and background
    n_threshold: proportion of background sequences that can be N
    genome: genome from which to draw background sequences
    """
    chromosomes = sorted(chrom_seq_dict.keys())
    numChromosomes = len(chrom_seq_dict.keys()) # number of chromosomes

    ### calculate GC content and average length of the target sequences
    target_gc_count = 0
    target_length_count = 0
    for pos in target_positions:
        seq = pos[-2]
        if len(seq) >0:
            target_gc_count += seq.count('G')
            target_gc_count += seq.count('C')
            target_length_count += len(seq)
    target_gc_content = (target_gc_count + 0.1)/(target_length_count+0.1) # GC content of target sequences
    
    ### select random genomic loci such that they do no overlap target sequences
    numSelected = 0
    # candidate pool of background seqs is size_ratio X larger
    numToSelect = len(target_positions) * size_ratio 
    candidate_positions = []
    numNallowed = int(n_threshold * mean_target_length) # number of allowable Ns
    counter = 0
    while numSelected < numToSelect:
        if counter % 100000 == 0:
            print(counter, numSelected)
        # select random chromsome
        chromIndex = np.random.randint(numChromosomes)
        randChrom = chromosomes[chromIndex]
        randChromSize = chrom_size_dict[randChrom]
        # must find non overlapping segment on this chromosome before moving on
        selectedSequence = False
        while not selectedSequence:
            counter += 1
            randStart = np.random.randint(randChromSize)
            randEnd = randStart + mean_target_length
            overlap_sum = np.sum(target_chr_position_dict[randChrom][randStart:(randEnd)])
            
            if not overlap_sum > 0:
                randSeq = chrom_seq_dict[randChrom][randStart:(randEnd+1)]
                numN = randSeq.count('N')
                if numN <= numNallowed:
                    rand_gc_count = randSeq.count('G')+ randSeq.count('C')
                    rand_gc = rand_gc_count/mean_target_length
                    if abs(target_gc_content - rand_gc) <= tolerance:
                        selectedSequence = True
                        numSelected+=1
                        candidate_positions.append([randChrom, randStart, randEnd, randSeq])
            if counter > 10000:
                break
    # calcuate GC content of background samples
    background_gc_count = 0
    background_length = 0
    for cp in candidate_positions:
        s = cp[3]
        background_gc_count += s.count('G')
        background_gc_count += s.count('C')
        background_length += len(s)
    background_gc_content = background_gc_count/(background_length+0.0000001)
    print('target GC:', target_gc_content, 
          'background GC:', background_gc_content, 
          'target length:', mean_target_length,
          'numTargetPositions',len(target_positions),
          'backgroundPositions', len(candidate_positions))
    return candidate_positions

def write_background_positions(background_positions, output_dir):
    """
    converts background positions into a bed file and a 
    """

    bed_file = open(output_dir + '/background.bed', 'w')
    fasta_file = open(output_dir + '/background.fasta', 'w')
    counter = 0
    for pos in background_positions:
        chrom = pos[0]
        start = str(pos[1])
        end = str(pos[2])
        seq = str(pos[3])
        randID = 'bg_' + str(np.random.randint(100000)) + '_' + str(counter)
        counter += 1
        bed_file.write('\t'.join([chrom, start, end, randID, '\n']))
        fasta_file.write('>' + randID + '\n')
        fasta_file.write(seq + '\n')
    bed_file.close()
    fasta_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Constructs random GC matched '+
                                     'background regions')
    parser.add_argument("inputPath",
        help="path to a bed file containing a chr, start, end, and strand column",
        type = str)
    parser.add_argument("genome",
        help="genome from which to construct background regions",
        type=str)
    parser.add_argument("outputPath",
        help="directory where output files should be written",
        default="./", type=str)
    parser.add_argument("-sizeRatio",
        help="size of the background region with respect to the target region",
        default = 1.0, type=float)
    parser.add_argument("-numBins",
        help="number of bins to use for GC normalization",
        default = 10,
        type=float)
    parser.add_argument("-nThreshold",
        help="maximum fraction of background sequences that can be N",
        default = 0.1,
        type=float)
    parser.add_argument("-filterChromosomes",
        help="chromosomes to ignore",
        type=str,
        default=['chrM', 'chrY'],
        nargs='+')

    # parse arguments
    args = parser.parse_args()

    input_path = args.inputPath
    output_path = args.outputPath
    size_ratio = args.sizeRatio
    num_bins = args.numBins
    n_threshold = args.nThreshold
    genome = args.genome
    filter_chromosomes = args.filterChromosomes
    
    target_positions = read_target_positions(input_path, filter_chromosomes)
    
    background_positions = get_random_background(target_positions, 
                                                size_ratio = size_ratio, 
                                                num_bins = num_bins, 
                                                n_threshold = n_threshold,
                                                genome = genome,
                                                filter_chromosomes=filter_chromosomes
                                                )
    write_background_positions(background_positions, output_path) 

