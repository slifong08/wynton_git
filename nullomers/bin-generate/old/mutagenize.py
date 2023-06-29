"""
Take kmer spectra, nullomer spectra, mutagenize sequences to get kmer spectra of the mutants. 

input
    mer (str) - identity of the kmer length
    nmuts (str) - number of positions to mutagenize simultaneously. 
    config (str) - path to config.ini file (must be in current working directory)
    
"""
from Bio.SeqUtils import gc_fraction

import gzip
from itertools import product, combinations
import numpy as np
import os, sys
import subprocess as sp
sys.path.append("/wynton/home/ahituv/fongsl/tools/py_")
import config_readwrite as crw

MER = sys.argv[1]  #str - 11, 12, 13 14mer etc 
NMUTS= int(sys.argv[2])  # number of mutations to test for each nullomer. 
config_tag = sys.argv[3]  #config_tag = "config"

config, cfn = crw.read_config(os.path.join(os.getcwd(), config_tag))
NULLS = config[f"{MER}mer"]["nullomers"]
KMERS = config[f"{MER}mer"]["kmers"]
RE= config[f"{MER}mer"]["results"]
PATH=config[f"{MER}mer"]["path"]

MER, NMUTS, = int(MER), int(NMUTS) # make string into int


###
# FUNCTIONS
###
def getSeqDict(file):
    """
    get dictionary of kmer sequences and counts, guides
    
    input
        file (str) - file w/ full path
    
    method 
        1. instantiate collection dict
        2. read file, extract sequence and kmer counts
    return
        seq_dict (dict) - dictionary of kmers and counts

    """
    print("getting sequence file", file)
    
    #1
    seq_dict = {}
    
    #2
    with gzip.open(file, "rt") as reader:
        for line in reader.readlines():
            seq, count = line.strip("\n").split(",")
            seq_dict[seq]=count
            
    return seq_dict

def genKmers(length):
    """
    return all sequence permutations, including repeats (AAAAA, GGGGG, CCCCC etc.)
    
    require
        itertools.product
    
    input
        length (int)
    
    return
        mers_list (list) - list of all nucleotide permutations 
    """
    print("generating kmer space length", length)
    
    mers = product("ACTG", repeat=length)
    
    mers_list = list("".join(i) for i in mers)
    
    return mers_list


def getPosLetterCombos(nmuts, sequence):
    """
    return combinations of (1) indices (2) mutated bases for mutating a sequence
    
    require 
        itertools
        
    inputs 
        nmuts (int) - number mutations to make
        sequence (str) - sequence to be mutated
        
    method
        1. get index combinations based on sequence length and number of mutations to make. This makes a map of all possible combinations of sequences to mutate
            1.1 combinations requires that each index is unique and non-redundant. 
                Order does not matter - 
                    e.g. (2,4,5) is the same as (5,2,4) because indexes 2, 4, and 5, will all be mutated.
                    
        2. get sequence product to mutate at indexes
            2.1 - product allows for repeats of the same base in different positions
        
    return
        mut_pos (list) - list of positional combinations to mutate
        mut_bases (list) - list of letter combinations to mutate
        
    
    """
    #print("making index combinations, nucletide permutations of length", nmuts)
    
    #1 index combinations
    mut_pos = list(combinations(np.arange(len(sequence)), nmuts))
    
    #2 nucleotide permutations per index combo. 
    mut_bases = list(product("ACGT", repeat=nmuts))
    
    return mut_pos, mut_bases

def buildSeqMut(sequence, mut_pos, mut_bases):
    
    """
    mutate sequence at position with letter
    multiple positions and letters can be inserted into the sequence simultaneously. 

    input
        sequence (str) - original sequence
        mut_pos (set) - Sets of single, tuple, threeple positional index(es) to mutate
        mut_bases (tuple) - Sets of single, tuple, threeple nucleotide combinations to mutate sequence to

    method
        1. instantiate seqs set to collect mutated sequences, add identity to seq set
        2. per positions to mutate in set
            2.1 per base combinations to mutate at these positions
        3. zip and iterate through positions and bases, mutating input sequence
        4. IF mut_seq != input sequence, then return. Else, skip

            
    return 
        seqs (set) - set of sequences with mutations 
    
    """

    
    #1
    seqs, mut_seq = set(), ""
    
    
    #2
    for pos in mut_pos:
        
        #2.1
        for letters in mut_bases:
            
            #3
            for p, l in zip(pos, letters):

                if mut_seq =="":
                    mut_seq = sequence[:p] + l + sequence[p + 1:]
                    
                else:
                    mut_seq = mut_seq[:p] + l + mut_seq[p + 1:]
        
            #4
            if mut_seq != sequence:  
                seqs.add(mut_seq)
                mut_seq = ""


            else:
                #print('no mut', sequence, mut_seq
                mut_seq = ""
                pass
            
    seqs.add(sequence)
    return seqs

def prettifySeq(original, mut):
    
    """
    prettify sequence. 
        All matching bases will be written as "."
        All non matching bases will be written w mutated base identity. 
    """
    prettyseq=""
    for o, m in zip(original, mut):
        if o!=m: 
            prettyseq+=m
        else:
            prettyseq += "."
            
    return prettyseq

def generateMismatchSpectra(input_dict, nmuts, null_spectra, kmer_spectra):
    
    """
    input 
        input_dict (dict) - dictionary where keys are sequence strings and values are kmer 
        nmuts (int) - max number of mismatches to mutate each kmer sequence by
        null_spectra (dict) - dictionary of nullomer keys and their frequency count (0)
        kmer_spectra (dict) - dictionary of kmer keys and their frequency count (value)

    require
        getPosLetterCombos function
        buildSeqMut function
        prettifySeq function
        
    method
        1. instantiate collection dictionary.
        
        2. per nullomer seq in dictionary

        3. get all combinations of indexes and mutated bases to try. 
            mutate bases will be a single, tuple, or threeple depending on the number of mutations desired.
            
        4. mutate all positions of the nullomer with that base. 

            4.1 - get mutated sequences using buildSeqMut function
                NOTE - mutBase function will remove identity sequence. See function method step 4.3.
        5. instantiate kmer_mismatch dict. Collect kmer distribution for each base mismatch. 
        6. per mutated sequence 
            6.1 look up kmer count of the mutated sequence. Test whether it is a prime. 
            6.2 add any resurfaced nullomer sequences to the resurface dict 
                - that is, resurfaced nullomers are sequences where mutations create nullomers
                - return this dictionary
            6.3 add kmer count to the dictionary as pretty key
        7. Add seq kmer dictionary back into the collection dictionary. 

    
    return
        collection_dict (dictionary) - summarizes every sequence in input dictionary with 
            every mismatch combination for each sequence 
            and per mismatch combo, the kmer counts for each mismatch

    """
    print("making kmer mismatch spectra w/ N mutatiosn =", nmuts)
    #1
    collection_dict, resurface_dict = {}, {}

    #2
    for i, seq in enumerate(input_dict.keys()):
        
        gc_seq = gc_fraction(seq)

        collection_dict[seq] = {}  # collect results for all mismatches per sequence. 

        #3
        mut_pos, mut_bases = getPosLetterCombos(nmuts, seq)

        #4
        seqs = buildSeqMut(seq, mut_pos, mut_bases)

        #5
        kmer_mismatch=collection_dict[seq]  # get sequence-specific dictionary

        #6
        for mut_seq in seqs:

            # 6.1 - kmer frequency
            if mut_seq not in null_spectra.keys():
                kmer_count=int(kmer_spectra[mut_seq])
                
            # 6.2 - resurfaced nullomers 
            else:  # if not identity sequence + not in nullomer dictionary. 
                
                resurface_dict[seq] = []  # resurface sequence key will have a list of resurfaced nullomers
                resurf_list = resurface_dict[seq]  # get list
                resurf_list.append(mut_seq)  # append to list
                resurface_dict[seq] = resurf_list  # update dictionary
                
                #print('\n', i, '\nprime', s, "\nnull ", seq, "\n", pretty, "\n")
                kmer_count=0

            #6.3 - prettify and record kmer count
            
            pretty = prettifySeq(seq, mut_seq)  # pretty key

            kmer_mismatch[pretty] = kmer_count

        #7 
        collection_dict[seq] = kmer_mismatch  # update sequence specific dictionary each time. 

    return collection_dict, resurface_dict

def writeDict(dict_, out):
    """
    write dictionary as gzip file
    
    input
        dict_ (dict) - dictionary of nullomer sequence keys + list of associated nullomers
        out (str) - path to file to write (gzipepd)
    
    method
        1. open out file
        2. unpack dictionary
        3. count number of mismatch nullomers corresponding to this nullomer
            3.1 if count is equal to the number of nullomers required to make a prime, set prime variable to True
            
        4. write all values
        5. close writer

    """
    #1
    with gzip.open(out, "wt") as writer:
        #2
        for key, value in dict_.items():
            
            #3
            n_nullomers = len(value)
            
            #3.1
            prime=False
            if n_nullomers == len(key)*3:
                prime =True
                
            #4
            writer.write(f"{key}\t{','.join(value)}\t{n_nullomers}\t{prime}\n")
        #5
        writer.close()
        
    print("wrote results", out)
    
    
def main(argv):
    
    # open up and get nullomer, kmer sequence dictionaries. 
    null_seqs, kmer_seqs = getSeqDict(NULLS), getSeqDict(KMERS)
    
    # generate mismatch spectra
    null, prime_dict = generateMismatchSpectra(null_seqs, NMUTS, null_seqs, kmer_seqs)
    
    # write dictionary
    OUT_NORDER = os.path.join(PATH, f"order.{NMUTS}.{MER}.nullomers.tsv.gz")
    writeDict(prime_dict, OUT_NORDER)
    
    config[f"{MER}mer"][f"{NMUTS}orderNull"] = OUT_NORDER
    crw.write_config(config, cfn)
    

    # get GC content to match kmer spectra on.
    gc = []

    for seq in null.keys():
        gc.append(gc_fraction(seq))

        
    # randomly sample 100x kmers. 
    random_kmer_dict = {}
    random_kmers = np.random.choice(list(kmer_seqs.keys()), size=(len(null_seqs.keys())*100))
    
    # get randomly sampled kmer spectra
    for i in random_kmers:
        if gc_fraction(i) in gc:
            random_kmer_dict[i] = 0

    kmer, prenullomer_dict = generateMismatchSpectra(random_kmer_dict, NMUTS, null_seqs, kmer_seqs)

    ### WRITE nullomer, pre-nullomer dictionary. 
    OUT_KORDER = os.path.join(PATH, f"order.{NMUTS}.{MER}.nullomers.tsv.gz")
    writeDict(prenullomer_dict, OUT_KORDER)
    
    config[f"{MER}mer"][f"{NMUTS}orderK"] = OUT_KORDER
    crw.write_config(config, cfn)
    
    
if __name__ =="__main__":
    main(sys.argv[1:])