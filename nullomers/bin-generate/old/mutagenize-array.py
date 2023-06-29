
"""
Take kmer spectra, nullomer spectra, mutagenize sequences to get kmer spectra of the mutants. 

input
    run (str) 
    mer (str) - identity of the kmer length
    nmuts (str) - number of positions to mutagenize simultaneously. 
    config (str) - path to config.ini file (must be in current working directory)
    
"""
import argparse
import os, sys
sys.path.append("/wynton/home/ahituv/fongsl/tools/py_")

import config_readwrite as crw
import glob
import gzip
from itertools import product, combinations

import numpy as np

import subprocess as sp

from timeit import default_timer as timer

parser = argparse.ArgumentParser()

parser.add_argument('arrayN', type=int, help="array number from job_id")
parser.add_argument('kmerLen', type=int, help="kmer length")
parser.add_argument('nMutations', type=int, help="n mutations to make to nullomer")
parser.add_argument('config', type=str, help = "config file for genome build, path to datasets")

args = parser.parse_args()
RUN = args.arrayN
MER = args.kmerLen
NMUTS= args.nMutations
config_tag = args.config

#RUN, MER, NMUTS, config_tag = "14", "1", "config.hg38"

config, cfn = crw.read_config(os.path.join(os.getcwd(), config_tag))
PATH = config[f"{MER}mer"]["path"]
ARRAY = config[f"{MER}mer"]["array"]
KMERS = config[f"{MER}mer"]["kmers"]
RE= config[f"{MER}mer"]["results"]
PATH=config[f"{MER}mer"]["path"]

# functions 

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
            
        4. write ONLY NULLOMER values
        5. close writer

    """
    #1
    with open(out, "wt") as writer:
        #2
        for key, value in dict_.items():

            #3
            n_nullomers = int(value)

            #3.1
            prime=False
            
            if n_nullomers == (len(key)*3)+1:  #3 for the other 3 bases + identity
                prime =True

                #4
                writer.write(f"{key}\t{n_nullomers}\t{prime}\n")

        #5
        writer.close()
        
    print("wrote results", out)

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

def generateMismatchSpectra(nullfile, nmuts, kmer_file):
    
    """
    input 
        null_file (str) - nullomer file w two columns (dictionary-like) where keys are sequence strings and values are kmer 
        nmuts (int) - max number of mismatches to mutate each kmer sequence by
        kmer_spectra (dict) - dictionary of kmer keys and their frequency count (value)

    require
        getPosLetterCombos function
        buildSeqMut function
        prettifySeq function
        
    method
        0. get name, nullomer dict, kmer dict
        
        1. instantiate resurfacing dictionary.
        
        2. per nullomer seq in dictionary (e.g. "AAA,0")

        3. get all combinations of indexes and mutated bases to try. 
            number of simultaneously mutated bases can be adjusted
            
        4. mutate all positions of the nullomer with that base.
            e.g. (AAA, AAC, AAG, AAT, AAA, ACA, AGA, ATA, AAA, CAA, GAA, TAA) 
            NOTE - identity (AAA) will be in the mutated product. 

            4.1 - get mutated sequences using buildSeqMut function
                NOTE - buildSeqMut function will skip identity sequence during processing, but add back in once instance of the identity sequence at the end. It is important to keep track of mutations that create the identity nullomer sequence.  See function method step 4.3.
        
        5. per mutated sequence 
            5.1 look up kmer count of the mutated sequence. If mut_seq is not kmer, assume mut_seq is nullomer, else continue
            5.2 add any resurfaced nullomer sequences to the resurface dict 
                - that is, resurfaced nullomers are sequences where mutations create nullomers
                - count the number of times a nullomer resurfaces. 
                - add nullomer resurfacing count to dictionary
        6. write dictionary to file. 

    
    return
        collection_dict (dictionary) - summarizes every sequence in input dictionary with 
            every mismatch combination for each sequence 
            and per mismatch combo, the kmer counts for each mismatch

    """
    #0
    print("null file", nullfile)
    
    name = (os.path.split(nullfile)[1]).strip(".csv.gz")
    null_seqs, kmer_spectra = getSeqDict(nullfile), getSeqDict(kmer_file)
    
    
    print("making kmer mismatch spectra w/ N mutations =", nmuts)
    #1
    resurface_dict = {}

    #2
    for i, seq in enumerate(null_seqs.keys()):
       
        #3
        mut_pos, mut_bases = getPosLetterCombos(nmuts, seq)

        #4
        seqs = buildSeqMut(seq, mut_pos, mut_bases)


        #5
        for mut_seq in seqs:

            # 5.1 - kmer frequency
            if mut_seq in kmer_spectra.keys():
                continue
                
            # 5.2 - nullomer resurfacing (i.e. when mutation produces yet another nullomer sequence)
            else:  # if not identity sequence + not in nullomer dictionary. 
                if seq not in resurface_dict:
                    resurface_dict[seq] = 1  # resurface sequence key will have a set of resurfaced nullomers
                else:
                    resurface_dict[seq] +=1
                
            
    #7 write dictionary to file
    out = os.path.join(PATH, f"order.{name}.tsv")
    
    print(len(resurface_dict.keys()), out)
    writeDict(resurface_dict, out)


def readArray(array, run_number):
    # get the array file to run. 
    with open(array, "r") as reader:
        for line in reader:
            n, file, key = (line.strip("\n")).split("\t")
            
            if int(n)==int(run_number):
                return file

# main

def main(argv):

    NULL = readArray(ARRAY, RUN)
    print("NULL file", NULL, "for run", RUN)

    # generate mismatch spectra, write n-order nullomer dictionary to file. 
    # e.g. for each file ALL.AAAA.14mers-nullomers.csv.gz, ALL.AAAC.14mers-nullomers.csv.gz, ALL.AAAG.14mers-nullomers.csv.gz
    generateMismatchSpectra(NULL, NMUTS, KMERS)
   
    

if __name__ =="__main__":
    main(sys.argv[1:])