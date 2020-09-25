import numpy as np
import pandas as pd

import gzip

from collections import defaultdict

from tqdm import tqdm
from os.path import join
from Bio import SeqIO


complementarity = {'A' : 'T',
                   'C' : 'G',
                   'G' : 'C',
                   'T' : 'A',
                   'N' : 'N'
                   }

def reverse_complement(text):
    '''
    This function returns the reverse complement ("reverse_text") of a DNA sequence ("text").
    '''
    N = len(text)
    reverse_text = ''
    for i in range(1, N+1):
        j = -1*i
        reverse_text += complementarity[text[j]]
    return reverse_text


def burrows_wheeler_transform(text):
    '''
    This function generates the BWT of text.
    '''
    bwt_matrix = [text[i:] + text[0:i] for i in range(len(text))]
    bwt_matrix.sort()
    bwt = ''
    for e in bwt_matrix:
        bwt += e[-1]
    return bwt


def partial_suffix_array_construction(text, K):
    '''
    This function creates a partial suffix array of text.
    '''
    suffixes = [(text[i:], i) for i in range(len(text))]
    suffixes.sort()
    partial_suffix_array = dict()
    for i in range(len(suffixes)):
        if suffixes[i][1]%K == 0:
            partial_suffix_array[i] = suffixes[i][1]
    return partial_suffix_array


def get_checkpoints(bwt, K):
    '''
    This function returns a checkpoint matrix from the Burrows-Wheeler Transform.
    '''
    counts = dict()
    unique_symbols = list(set(bwt))

    # Initialize Counts
    for symbol in unique_symbols:
        counts[(0, symbol)] = 0

    for i in range(1, len(bwt) + 1):
        for symbol in unique_symbols:
            counts[(i, symbol)] = 0

    # Get Counts
    for i in range(len(bwt)):
        next_symbol = bwt[i]
        for symbol in unique_symbols:
            if symbol == next_symbol:
                counts[(i + 1, symbol)] = counts[(i, symbol)] + 1
            else:
                counts[(i + 1, symbol)] = counts[(i, symbol)]

    # Get checkpoints
    checkpoints = dict()
    for key in counts.keys():
        if key[0] % K == 0:
            checkpoints[key] = counts[key]
    return checkpoints


def get_first_occurence(first_column_list):
    '''
    This function returns an array with the first occurrences of unique symbols in first column.
    '''
    first_occurrence = dict()
    for symbol in set(first_column_list):
        first_occurrence[symbol] = first_column_list.index(symbol)
    return first_occurrence


def better_bw_matching_checkpoints(first_occurrence, last_column, pattern, checkpoints, K, partial_suffix_array):
    '''
    This function performs the Burrows-Wheeler Matching of a pattern using checkpoints.
    '''
    top = 0
    bottom = len(last_column) - 1
    tmp_pattern = pattern
    while top <= bottom:
        if len(tmp_pattern) > 0:
            symbol = tmp_pattern[-1]
            tmp_pattern = tmp_pattern[:-1]

            # Update top
            if top%K == 0:
                top = first_occurrence[symbol] + checkpoints[(top, symbol)]
            else:
                pos = (top // K) * K
                new_top = first_occurrence[symbol] + checkpoints[(pos, symbol)]
                for i in range(pos, top):
                    if last_column[i] == symbol:
                        new_top += 1
                top = new_top

            # Update bottom
            if (bottom + 1)%K == 0:
                bottom = first_occurrence[symbol] + checkpoints[(bottom + 1, symbol)] - 1
            else:
                pos = ((bottom + 1) // K) * K
                new_bottom = first_occurrence[symbol] + checkpoints[(pos, symbol)]
                for i in range(pos, bottom + 1):
                    if last_column[i] == symbol:
                        new_bottom += 1
                bottom = new_bottom - 1
        else:
            break
    if top > bottom:
        return []
    else:
        # Look for starting positions using partial suffix array
        positions = []
        for i in range(top, bottom + 1):
            if i in partial_suffix_array.keys():
                positions.append(partial_suffix_array[i])
            else:
                top = i
                extra_steps = 0
                while top not in partial_suffix_array.keys():
                    symbol = last_column[top]

                    # Update top
                    if top % K == 0:
                        top = first_occurrence[symbol] + checkpoints[(top, symbol)]
                    else:
                        pos = (top // K) * K
                        new_top = first_occurrence[symbol] + checkpoints[(pos, symbol)]
                        for i in range(pos, top):
                            if last_column[i] == symbol:
                                new_top += 1
                        top = new_top
                    extra_steps += 1
                positions.append(partial_suffix_array[top] + extra_steps)
    
    # Add 1 to each position
    positions = [p+1 for p in positions]
    return positions


def multiple_better_bwmatching_checkpoints(patterns, text, K):
    '''
    This function performs the bw_matching for a set of patterns. Then, it returns a list of number of matches for all
    of them. *** In this case, partial suffix array and checkpoints are used. ***
    Based on ROSALIND format of inputs.
    '''
    # Symbols
    text_symbols = set(text)
    
    # Add $
    text = text + '$'

    # BWT
    bwt = burrows_wheeler_transform(text)

    # First Occurrence
    first_column_list = [char for char in bwt]
    first_column_list.sort()

    first_occurrence = get_first_occurence(first_column_list)
    del first_column_list

    # Checkpoints
    checkpoints = get_checkpoints(bwt, K)

    # Partial suffix array
    partial_suffix_array = partial_suffix_array_construction(text, K)

    # Matches
    matches = dict()
    
    for pattern in patterns:
        pattern_symbols = set(pattern)
        if (pattern_symbols&text_symbols) == pattern_symbols:
            aligns = better_bw_matching_checkpoints(first_occurrence, bwt, pattern, checkpoints, K, partial_suffix_array)
            matches[pattern] = set(aligns)
        else:
            matches[pattern] = set()
    return matches


def search_positions(patterns, ref_seq, K=100):
    '''
    This function performs the BWA of a set of patterns on a reference sequence, considering its 
    complementary sequence.
    '''
    # Perform search in ref_seq
    ref_align = multiple_better_bwmatching_checkpoints(patterns, ref_seq, K)
    
    # Perform search in complementary seq
    comp_seq = reverse_complement(ref_seq)
    comp_align = multiple_better_bwmatching_checkpoints(patterns, comp_seq, K)
    
    results = []
    # Positive strand
    for k, v in ref_align.items():
        for p in v:
            start = p
            end = start-1+len(k)
            record = (k, int(start), int(end), '+')
            results.append(record)
    
    # Negative strand
    N = len(ref_seq)
    for k, v in comp_align.items():
        if len(v) > 0:
            for p in v:
                start = N+1-(p+len(k))
                end = N+1-p
                record = (k, int(start), int(end), '-')
                results.append(record)
        else:
            record = (k, np.nan, np.nan, np.nan)
            results.append(record)
    df = pd.DataFrame.from_records(results, columns=['Pattern', 'Start-Position', 'End-Position', 'Strand'])
        
    return df


def read_gff(filename, file_format='gzip'):
    '''
    This function reads a GFF files.
    
    Parameters
    __________
    
    filename : str
        It is the path to the GFF file.
    file_format : str
        It could be gzip or uncompress. Here, we use gzip by default.
        
    Returns
    _______
    
    df : DataFrame
        Pandas dataframe containing gff info.
    '''
    entries = []
    if file_format == 'gzip':
        with gzip.open(filename, 'rb') as f:
            for i, line in enumerate(f):
                if i % 10000 == 0:
                    print("Reading line {}".format(i))
                    
                line = line.decode()
                #if line.startswith('#'):
                if line[0] == '#':
                    continue 

                vals = line.split('\t')
                seqname = vals[0]
                feature = vals[2]
                start = int(vals[3])
                end = int(vals[4])
                strand = vals[6]
                idx = vals[8].split(';').lstrip('ID=')
                
                
                start, end = int(start), int(end)
                
                entries.append((idx, feature, start, end, strand, seqname))
    else:
        with open(filename, 'rb') as f:
            for i, line in enumerate(f):
                if i % 10000 == 0:
                    print("Reading line {}".format(i))
                    
                #if line.startswith('#'):
                line = line.decode()
                if line[0] == '#':
                    continue 

                vals = line.split('\t')
                seqname = vals[0]
                feature = vals[2]
                start = int(vals[3])
                end = int(vals[4])
                strand = vals[6]
                idx = vals[8].split(';')[0].lstrip('ID=')
                
                
                start, end = int(start), int(end)
                
                entries.append((idx, feature, start, end, strand, seqname))
                
    df = pd.DataFrame.from_records(entries, columns=['Id', 'Type', 'Start', 'End', 'Strand', 'Location'])
    #df.set_index('Id', inplace=True)
    del entries
    return df


def map_element_with_gff(element_start, element_end, gff_df, strand=None):
    mapped_idxs = []
    for row in gff_df.iterrows():
        idx = row[0]
        r = row[1]
        if (int(r['Start']) <= element_start) & (int(r['End']) >= element_end):
            if strand is None:
                mapped_idxs.append(idx)
            else:
                if strand != r['Strand']:
                    mapped_idxs.append(idx)                  
            
    return gff_df.loc[mapped_idxs, :].reset_index(drop=True)


def between(a, b, x):
    if (a <= x) & (x <= b):
        return True
    else:
        return False
    
    
def interval_extract(list): 
    '''Taken from https://www.geeksforgeeks.org/python-make-a-list-of-intervals-with-sequential-numbers/'''
    length = len(list) 
    i = 0
    while (i< length): 
        low = list[i] 
        while i <length-1 and list[i]+1 == list[i + 1]: 
            i += 1
        high = list[i] 
        if (high - low >= 1): 
            yield [low, high] 
        elif (high - low == 1): 
            yield [low, ] 
            yield [high, ] 
        else: 
            yield [low, ] 
        i += 1
        
        
from Bio import SeqIO

def get_read_sequences(filename, file_format='fastq', compression=None):
    '''Open FASTQ file and get all sequences in that file.'''
    print('Loading {}'.format(filename))
    if compression == None:
        with open(filename, "rU") as handle:
            reads = dict()
            for record in SeqIO.parse(handle, file_format):
                reads[str(record.id)] = str(record.seq)
    elif compression == 'gzip':
        with gzip.open(filename, "rt") as handle:
            reads = dict()
            for record in SeqIO.parse(handle, file_format):
                reads[str(record.id)] = str(record.seq)
    return reads