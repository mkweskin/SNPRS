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

def windowsPool(pool_cpu,raw_alignment,raw_spp,bases,compare):
    pool = mp.Pool(pool_cpu)
    var_getPairwise=partial(getPairwise_Windows, raw_alignment = raw_alignment,raw_spp = raw_spp, bases = bases)
    pairwise_results = pool.map(var_getPairwise, compare_list)
    pool.close()
    return(pairwise_results)

def getPairwise_Windows(raw_alignment,raw_spp,bases,compare):
    # Populate temp alignment with two taxa from compare
    temp_alignment = raw_alignment[0:0]
    for i in compare:
        temp_alignment.add_sequence(str(raw_alignment[raw_spp.index(i)].id), str(raw_alignment[raw_spp.index(i)].seq))

    # If resulting alignment is empty, raise exception
    if int(temp_alignment.get_alignment_length()) == 0:
        sys.exit("ERROR: Alignment processed, but appears to have no bases...")
    else:
        align_length = temp_alignment.get_alignment_length()

    # Get columns where both species have bases, and if so, find identical sites
    cols_with_base = []
    identical_sites = []

    for i in range(0,align_length):
        seq_string = temp_alignment[:, i]
                
        # Get alleles
        allele_list = list(Counter(seq_string))
        allele_count = len(allele_list)
        
        if all(b in bases for b in seq_string):
            cols_with_base.append(i)
            
            if allele_count == 1:
                identical_sites.append(i)
            else:
                continue
        else:
            continue

    compare_id = ";".join(sorted([compare[0],compare[1]]))
    if len(cols_with_base) == 0:
        if not len(identical_sites) == 0:
            sys.exit("ERROR: Identical sites found but no co-called?")
        else:
            return(pd.DataFrame([[compare[0],compare[1],compare_id,0,0,0]],columns=['Sample_A','Sample_B','Comparison','Cocalled_Sites','Identical_Sites','Prop_Identical']))
    else:
        return(pd.DataFrame([[compare[0],compare[1],compare_id,len(cols_with_base),len(identical_sites),float(len(identical_sites))/float(len(cols_with_base))]],columns=['Sample_A','Sample_B','Comparison','Cocalled_Sites','Identical_Sites','Prop_Identical']))

def getPairwise(compare):
    global raw_alignment
    global bases
    global raw_spp
    
    # Populate temp alignment with two taxa from compare
    temp_alignment = raw_alignment[0:0]
    for i in compare:
        temp_alignment.add_sequence(str(raw_alignment[raw_spp.index(i)].id), str(raw_alignment[raw_spp.index(i)].seq))
    
    # If resulting alignment is empty, raise exception
    if int(temp_alignment.get_alignment_length()) == 0:
        sys.exit("ERROR: Alignment processed, but appears to have no bases...")
    else:
        align_length = temp_alignment.get_alignment_length()
    
    # Get columns where both species have bases, and if so, find identical sites
    cols_with_base = []
    identical_sites = []
    
    for i in range(0,align_length):
        seq_string = temp_alignment[:, i]
            
        # Get alleles
        allele_list = list(Counter(seq_string))
        allele_count = len(allele_list)
        
        if all(b in bases for b in seq_string):
            cols_with_base.append(i)
            
            if allele_count == 1:
                identical_sites.append(i)

    compare_id = ";".join(sorted([compare[0],compare[1]]))
    
    if len(cols_with_base) == 0:
        if not len(identical_sites) == 0:
            sys.exit("ERROR: Identical sites found but no co-called?")
        else:
            return(pd.DataFrame([[compare[0],compare[1],compare_id,0,0,0]],columns=['Sample_A','Sample_B','Comparison','Cocalled_Sites','Identical_Sites','Prop_Identical']))
    else:
        return(pd.DataFrame([[compare[0],compare[1],compare_id,len(cols_with_base),len(identical_sites),float(len(identical_sites))/float(len(cols_with_base))]],columns=['Sample_A','Sample_B','Comparison','Cocalled_Sites','Identical_Sites','Prop_Identical']))

def pairwiseAlign(path_to_align):
    
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
        
    # Create comparison list
    spp_count = len(raw_spp)
    compare_list = list()

    for sample_a in range(0,spp_count-1):
        for sample_b in range(sample_a+1,spp_count):
            compare_list.append((sorted([raw_spp[sample_a],raw_spp[sample_b]])))            

    # Detect CPUs and create a pool of workers
    if mp.cpu_count() == 1:
        pool_cpu = 1
    else:
        pool_cpu = mp.cpu_count() - 1    
    
    if os.name == 'nt':
        pairwise_results = windowsPool(pool_cpu,raw_alignment,raw_spp,bases,compare)
    else:
        with mp.Pool(pool_cpu) as pool:
            pairwise_results = pool.map(getPairwise, compare_list)
    
    pairwise_results = pd.concat(pairwise_results)
    pairwise_results.to_csv(alignment_path+'_Pairwise.tsv', sep = '\t', index=False,header=True)

def main(path_to_align):
    pairwiseAlign(path_to_align)

if __name__ == "__main__":
    main(sys.argv[1])
