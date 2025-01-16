#!/usr/bin/env python3
"""
    This script will take an alignment file (NEXUS, FASTA, or Phylip-relaxed) and:

    1) Prune down to a supplied subset of taxa
    2) Breakdown the signal at each site
"""

import sys
import os
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO
from operator import itemgetter
from collections import Counter
import multiprocessing as mp
from itertools import chain
import copy
import pickle
import tempfile
from functools import partial

def findOccurrences(s, ch):
    # findOccurrences returns indexes of occurrences of list items in a list or string
    return [i for i, letter in enumerate(s) if letter in ch]
    
def windowsPool(pool_cpu,bases,raw_alignment,raw_spp,reference_species,classify_samples):
    pool = mp.Pool(pool_cpu)
    var_getClassifier=partial(getClassifier_Windows, raw_alignment = raw_alignment,raw_spp = raw_spp, bases = bases,reference_species = reference_species)
    classifier_results = pool.map(var_getClassifier,classify_samples)
    pool.close()
    return(classifier_results)

def getClassifier_Windows(bases,raw_alignment,raw_spp,reference_species,classify):

    # Create a reference dictionary for matches to each reference species, as well as singletons
    conditions = reference_species + ['Singleton','Background']
    classify_dict = { i : 0 for i in conditions }
    classifier_sites = 0

    # Populate temp alignment with two taxa from compare
    temp_alignment = raw_alignment[0:0]
    all_spp = [classify] + reference_species
    for i in all_spp:
        temp_alignment.add_sequence(str(raw_alignment[raw_spp.index(i)].id), str(raw_alignment[raw_spp.index(i)].seq))

    # If resulting alignment is empty, raise exception
    if int(temp_alignment.get_alignment_length()) == 0:
        sys.exit("ERROR: Alignment processed, but appears to have no bases...")
    else:
        align_length = temp_alignment.get_alignment_length()

    # Go site by site, and identify positions where the test sample matches with a unique reference allele    
    for j in range(0,align_length):

        seq_string = temp_alignment[:, j]
        if all(b in bases for b in seq_string):
            sample_base = seq_string[0]
            sample_base_index = findOccurrences(seq_string, sample_base)
            
            # If the sample matches no reference samples, note as singleton
            if len(sample_base_index) == 1:
                classify_dict['Singleton'] = classify_dict.get("Singleton", 0) + 1
                
            # If the sample base matches with a unique reference base, update the dictionary
            if len(sample_base_index) == 2:
                reference_index = sample_base_index[1]
                reference_match = all_spp[reference_index]
                classify_dict[reference_match] = classify_dict.get(reference_match, 0) + 1
                classifier_sites += 1

            # If the sample base matches with more than one reference, update the dictionary
            if len(sample_base_index) > 2:
                classify_dict['Background'] = classify_dict.get("Background", 0) + 1

    # Convert dictionary to df
    sample_df = pd.DataFrame(list(classify_dict.items()),columns = ['Condition','Count'])
    reference_site_cols = ['Count_' + str(x) for x in sample_df['Condition']]
    sample_df['Condition'] = reference_site_cols
    sample_df = sample_df.transpose()
    sample_df.columns = sample_df.iloc[0]
    sample_df = sample_df.iloc[1:]
    sample_df['All_Sites'] = classifier_sites + classify_dict['Singleton'] + classify_dict['Background']
    sample_df['Classifier_Sites'] = classifier_sites
    sample_df['Sample_ID'] = classify
    
    return(sample_df)
        
def getClassifier(classify):
    global bases
    global raw_alignment
    global raw_spp
    global reference_species
    
    # Create a reference dictionary for matches to each reference species, as well as singletons
    conditions = reference_species + ['Singleton','Background']
    classify_dict = { i : 0 for i in conditions }
    classifier_sites = 0

    # Populate temp alignment with two taxa from compare
    temp_alignment = raw_alignment[0:0]
    all_spp = [classify] + reference_species
    for i in all_spp:
        temp_alignment.add_sequence(str(raw_alignment[raw_spp.index(i)].id), str(raw_alignment[raw_spp.index(i)].seq))

    # If resulting alignment is empty, raise exception
    if int(temp_alignment.get_alignment_length()) == 0:
        sys.exit("ERROR: Alignment processed, but appears to have no bases...")
    else:
        align_length = temp_alignment.get_alignment_length()

    # Go site by site, and identify positions where the test sample matches with a unique reference allele    
    for j in range(0,align_length):

        seq_string = temp_alignment[:, j]
        if all(b in bases for b in seq_string):
            sample_base = seq_string[0]
            sample_base_index = findOccurrences(seq_string, sample_base)
            
            # If the sample matches no reference samples, note as singleton
            if len(sample_base_index) == 1:
                classify_dict['Singleton'] = classify_dict.get("Singleton", 0) + 1
                
            # If the sample base matches with a unique reference base, update the dictionary
            if len(sample_base_index) == 2:
                reference_index = sample_base_index[1]
                reference_match = all_spp[reference_index]
                classify_dict[reference_match] = classify_dict.get(reference_match, 0) + 1
                classifier_sites += 1

            # If the sample base matches with more than one reference, update the dictionary
            if len(sample_base_index) > 2:
                classify_dict['Background'] = classify_dict.get("Background", 0) + 1

    # Convert dictionary to df
    sample_df = pd.DataFrame(list(classify_dict.items()),columns = ['Condition','Count'])
    reference_site_cols = ['Count_' + str(x) for x in sample_df['Condition']]
    sample_df['Condition'] = reference_site_cols
    sample_df = sample_df.transpose()
    sample_df.columns = sample_df.iloc[0]
    sample_df = sample_df.iloc[1:]
    sample_df['All_Sites'] = classifier_sites + classify_dict['Singleton'] + classify_dict['Background']
    sample_df['Classifier_Sites'] = classifier_sites
    sample_df['Sample_ID'] = classify
    
    return(sample_df)
        
def classifySamples(path_to_align,ref_species):

    # Set path to alignment
    alignment_path = str(path_to_align)

    global bases
    bases = ['A', 'C', 'T', 'G', 'a', 'g', 't', 'c','-']

    # Read in alignment and prune to desired species if requested
    global raw_alignment
    try:
        formats = {'nex': 'nexus', 'nexus': 'nexus',
                   'phy': 'phylip', 'phylip-relaxed': 'phylip-relaxed', 'phylip': 'phylip',
                   'fa': 'fasta', 'fasta': 'fasta'}

        fformat = formats[alignment_path.split('.')[-1]]
        raw_alignment = AlignIO.read(alignment_path, fformat)

    # If alignment cannot be read in, raise exception
    except:
        sys.exit("ERROR: Cannot process "+os.path.basename(alignment_path))

    # Get species from raw alignment
    global raw_spp
    raw_spp = list()
    for seq_record in raw_alignment:
        raw_spp.append(str(seq_record.id))
        
    # Check that all reference species are in the alignment
    global reference_species
    reference_species = sorted(str(ref_species).split(";"))    
        
    if not all(sp in raw_spp for sp in reference_species):
        sys.exit("ERROR: Specified reference species are not found in "+os.path.basename(alignment_path))
            
    # Get sample IDs to classify
    global classify_samples
    classify_samples = [x for x in raw_spp if x not in reference_species]

    # Detect CPUs and create a pool of workers
    if mp.cpu_count() == 1:
        pool_cpu = 1
    else:
        pool_cpu = mp.cpu_count() - 1

    if os.name == 'nt':
        classifier_results = windowsPool(pool_cpu,bases,raw_alignment,raw_spp,reference_species,classify_samples)
    else:
        with mp.Pool(pool_cpu) as pool:
            classifier_results = pool.map(getClassifier, classify_samples)

    classifier_results = pd.concat(classifier_results)
    reference_site_cols = ['Count_' + str(x) for x in reference_species]
    final_df_cols = ['Sample_ID'] + reference_site_cols + ['Classifier_Sites','All_Sites','Count_Singleton','Count_Background']

    classifier_results = classifier_results[final_df_cols]
    classifier_results.to_csv(alignment_path.replace("."+fformat,'_Classifier.tsv'), sep = '\t', index=False,header=True)

def main(path_to_align,ref_species):
    classifySamples(path_to_align,ref_species)

if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2])

