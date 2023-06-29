"""
Take kmer spectra, nullomer spectra, mutagenize sequences to get kmer spectra of the mutants. 

input
    jobnumber (int) - This number corresponds to the array element to be run. 
    array (str) - abs path to array file. 
   
"""
from timeit import default_timer as timer
import subprocess as sp
import numpy as np
from itertools import product, combinations
import gzip
import glob
from collections import Counter
import argparse
import os
import sys

parser = argparse.ArgumentParser()

parser.add_argument('narray', type=int, help="array number from job_id")
parser.add_argument('array', type=str, help="Array file")



args = parser.parse_args()
JOBNUM = args.narray
ARRAY = args.array


"""
JOBNUM = 1
ARRAY = "/wynton/home/ahituv/fongsl/nullomers/bin-generate/arrays/array-hs1.14mers.4.tsv"

"""


print("RUN", JOBNUM)


# functions

def array_reader(array, job_number):
    """
    read array, return values matching job number
    """

    with open(array, "r") as reader:
        for line in reader:
            num, key_pool, path, build, windowsize, keysize, writenull, nmuts = line.strip(
                "\n").split("\t")
            
            if int(num) == int(job_number):

                # change datatypes
                windowsize, keysize, writenull, nmuts = int(
                    windowsize), int(keysize), bool(writenull), int(nmuts)
                key_pool = key_pool.split(",")  # make list

                return key_pool, path, build, windowsize, keysize, writenull, nmuts


def getSeqDict(file_list):
    """
    get dictionary of kmer sequences and counts, guides

    input
        file_list (str) - files w/ full path

    method 
        1. instantiate collection dict
        2. read file, extract sequence and kmer counts
    return
        seq_dict (dict) - dictionary of kmers and counts
    """

    # 1
    seq_dict = {}

    if type(file_list) is str:  # if only one file, make it into a list w one item
        file_list = [file_list]

    print("getting sequence file", type(file_list))

    # 2
    for file in file_list:
        with gzip.open(file, "rt") as reader:
            for line in reader.readlines():
                if "nullomer" in file:
                    seq, count = line.strip("\n"), 0
                else:
                    seq, count = line.strip("\n").split(",")
                seq_dict[seq] = count

    return seq_dict


def genKmers(length):
    """
    return all sequence permutations, including repeats (AAAAA, GGGGG, CCCCC etc.)

    require
        itertools.product

    input
        length (int)
        
    method
        1. make product of ACTG sequences of length
        2. turn product into list

    return
        mers_list (list) - list of all nucleotide permutations 
    """
    
    print("generating possible kmer space length", length)
    
    #1
    mers = product("ACTG", repeat=length)

    #2
    mers_list = list("".join(i) for i in mers)

    return mers_list


def getPosLetterCombos(nmuts, kmer_len):
    """
    return combinations of (1) indices (2) mutated bases for mutating a sequence

    require 
        itertools

    inputs 
        nmuts (int) - number mutations to make
        kmer_len (int) - length of kmer sequence

    method
        1. Get index combinations based on sequence length and number of mutations to make. 
           This makes a map of all possible combinations of sequences to mutate
           
            1.1 combinations requires that each index is unique and non-redundant. 
                Order does not matter - 
                    e.g. (2,4,5) is the same as (5,2,4) because indexes 2, 4, and 5, will all be mutated.

        2. Get sequence product to mutate at indexes
            2.1 - product allows for repeats of the same base in different positions

    return
        mut_pos (list) - list of positional combinations to mutate
        mut_bases (list) - list of letter combinations to mutate


    """
    # print("making index combinations, nucletide permutations of length", nmuts)

    # 1 index combinations
    mut_pos = list(combinations(np.arange(kmer_len), nmuts))

    # 2 nucleotide permutations per index combo.
    mut_bases = list(product("ACGT", repeat=nmuts))

    return mut_pos, mut_bases


def buildSeqMut(sequence, mut_pos, mut_bases):
    """
    mutate sequence at position with mutant base
    multiple positions and bases can be inserted into the sequence simultaneously. 

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
        5. Add back in identity sequence


    return 
        seqs (set) - set of sequences with mutations 

    """

    # 1
    seqs, mut_seq = set(), ""

    # 2
    for pos in mut_pos:

        # 2.1
        for letters in mut_bases:

            # 3
            for p, l in zip(pos, letters):

                if mut_seq == "":
                    mut_seq = sequence[:p] + l + sequence[p + 1:]

                else:
                    mut_seq = mut_seq[:p] + l + mut_seq[p + 1:]

            # 4
            if mut_seq != sequence:
                seqs.add(mut_seq)
                mut_seq = ""

            else:
                # print('no mut', sequence, mut_seq
                mut_seq = ""
                pass

    # 5
    seqs.add(sequence)

    return seqs


def prettifySeq(original, mut):
    """
    prettify sequence. 
        All matching bases will be written as "."
        All non matching bases will be written w mutated base identity. 
    """
    prettyseq = ""
    for o, m in zip(original, mut):
        if o != m:
            prettyseq += m
        else:
            prettyseq += "."

    return prettyseq


def generateMismatchSpectra(nullfile, nmuts, kmer_spectra, path, windowsize):
    """
    input 
        null_file (str) - nullomer file w two columns (dictionary-like) 
                          where keys are sequence strings and values are kmer 
                          
        nmuts (int) - max number of mismatches to mutate each kmer sequence by
        kmer_spectra (dict) - dictionary of kmer keys and their frequency count (value)
        windowsize (int) - size of kmer
        path (str) - path to outdir

    require
        getPosLetterCombos function
        buildSeqMut function
        prettifySeq function

    method
        0. get name, nullomer dict

        1. instantiate resurfacing dictionary as Counter

        2. get all combinations of indexes and mutated bases to try. 
            number of simultaneously mutated bases can be adjusted
            - only needs to be made once for all sequences. 

        3. per nullomer seq in dictionary (e.g. "AAA,0")

        4. mutate all positions of the nullomer with that base.
            e.g. (AAA, AAC, AAG, AAT, AAA, ACA, AGA, ATA, AAA, CAA, GAA, TAA) 
            NOTE - identity (AAA) will be in the mutated product. 

            4.1 - get mutated sequences using buildSeqMut function
                NOTE - buildSeqMut function will skip identity sequence during processing, 
                        but add back in once instance of the identity sequence at the end. 
                     - It is important to keep track of mutations that create the identity nullomer sequence.  
                     - See buildSeqMut method step 5.

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
    # 0
    print("null file", nullfile)

    name = (os.path.split(nullfile)[1]).strip(".csv.gz")
    null_seqs = getSeqDict(nullfile)

    print("making kmer mismatch spectra w/ N mutations =", nmuts)

    # 1
    resurface_dict = Counter()

    # 2
    # set of mutation positions and tuples of mutated bases only needs to be made once
    mut_pos, mut_bases = getPosLetterCombos(nmuts, windowsize)

    # 3
    for i, seq in enumerate(null_seqs.keys()):

        # 4
        seqs = buildSeqMut(seq, mut_pos, mut_bases)

        # 5
        for mut_seq in seqs:

            # 5.1 - is nullomer kmer? 
            if mut_seq in kmer_spectra.keys():
                continue

            # 5.2 - nullomer resurfacing (i.e. when mutation produces yet another nullomer sequence)
            else:  # if not identity sequence + not in nullomer dictionary.
                resurface_dict[seq] += 1

    # 7 write dictionary to file
    out = os.path.join(path, f"order.{str(nmuts)}.{name}.tsv")

    print(len(resurface_dict.keys()), out)
    
    # function to write dictionary to the outfile
    writeDict(resurface_dict, out)


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
        4. write ONLY NULLOMERs
        5. close writer

    """
    # 1
    with open(out, "wt") as writer:
        
        # 2
        for key, value in dict_.items():

            # 3
            n_nullomers = int(value)

            # 3.1
            prime = False

            if n_nullomers == (len(key)*3)+1:  # 3 for the other 3 bases + identity
                prime = True

                # 4
                writer.write(f"{key}\n")

        # 5
        writer.close()

    print("wrote results", out)

    
# main

def main(argv):
    """
    Find n-order nullomers
    
    Take kmer spectra, nullomer spectra, mutagenize sequences to get kmer spectra of the mutants. 
    
    method
        1. read array and extract key variables
        2. for each prefix key in a pool of keys
        3. get the nullomer sequences, get the kmer sequences. 
            3.1 double check that none of the kmer files are nullomer files. 
        4. mutate nullomers NMUT times 
             determine if they are still nullomers after mutating NMUT times
             determine whether every NMUT number of mutations mutation creates a nullomer at every position. 
             If mutation in every position creates another nullomer sequence, this is also a nullomer. 
    """

    # read array. Get run info
    KEY_POOL, PATH, BUILD, WINDOW_SIZE, KEYSIZE, WRITENULL, NMUTS = array_reader(
        ARRAY, JOBNUM)

    for KEY in KEY_POOL:

        NULL = os.path.join(PATH, f"{KEY}.nullomers.csv.gz")

        KMER = glob.glob(os.path.join(PATH, f"*.csv.gz"))

        KMERS = []
        
        # filter out nullomer.csv files
        [KMERS.append(k) if "nullomer" not in k else next for k in KMER]

        # open the KMER dictionary once.
        KMER_DICT = getSeqDict(KMERS)

        # generate mismatch spectra, write n-order nullomer dictionary to file.
        # e.g. for each file ALL.AAAA.14mers-nullomers.csv.gz, 
                # ALL.AAAC.14mers-nullomers.csv.gz, ALL.AAAG.14mers-nullomers.csv.gz

        generateMismatchSpectra(NULL, NMUTS, KMER_DICT, PATH, WINDOW_SIZE)


if __name__ == "__main__":
    main(sys.argv[1:])