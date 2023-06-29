#!/usr/bin/env python
# coding: utf-8

"""
build prefix key files of specified key size. 
    scan genome for kmer with specified windowsize, step size=1
    count and write kmers for the prefix size. 
    
    if windowsize is small enough, generate potential kmer universe
    and explicitely write nullomer universe.

potential issue - need to download and update path to fa dictionary in extractFaSeq function. 
"""

import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from collections import Counter
from functools import partial
import glob
import gzip
from itertools import product
import numpy as np
import os
import sys

from timeit import default_timer as timer


arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("job_num", type=str, help='job number to look up in array.tsv')
arg_parser.add_argument("array", type=str, help='ARRAY file of job numbers, chromosome ids')

args = arg_parser.parse_args()
JOBNUM, ARRAY, =  args.job_num, args.array

"""
JOBNUM, ARRAY = 3, "/wynton/home/ahituv/fongsl/nullomers/bin-generate/arrays/array-sacCer3.11mer.tsv"
# PATH, WINDOW_SIZE, KEYSIZE, WRITENULL = "/wynton/home/ahituv/fongsl/dna/sacCer3", 11, 4, True
"""

# FUNCTIONS

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
                windowsize, keysize, writenull = int(
                    windowsize), int(keysize), bool(writenull)
                
                key_pool = key_pool.split(",")  # make list

                return key_pool, path, build, windowsize, keysize, writenull

# make outdirectory for results


def makeOutDir(path):
    """
    make outdirectory
    """

    print("OUTDIR", path)
    
    # if the outdir does not exist
    if os.path.exists(path) is False:
        
        # if the path to the outdir does not exist
        if os.path.exists(os.path.dirname(path)) is False:

            # make this directory first
            os.mkdir(os.path.dirname(path))

        # then make this directory second
        os.mkdir(path)
        
    return path


def makeKeys(keysize):
    """
    input
        keysize (int) - length of kmer-keys to create

    require
        itertools.product
        python list comprehension

    method
        itertools.product kmers -> join kmers into str -> make list

    return
        key_set (set) - list of strings. 
    """
    
    key_set = set()
    for item in product("ACGT", repeat=keysize):
        key_set.add("".join(item))

    return key_set


# ## extract fa sequence
def extractFaSeq(path, key_pool, keysize, window_size, build):
    """
    get chr-specific sequence, reverse complement

    require
        Bio.SeqIO.FastaIO - SimpleFastaParser
        gzip

    input
        chr_num (str) - chromosome number e.g. 1,2,3,X,Y
        path (str) - path to .fa files

    method
        1. get fa.gz name
        2. open fa file
        3. get sequence w/ simpleFastaParser
        4. get reverse complement by turning seq into Seq object

    return 
        seq (seq record) - str of sequence (should be one per chromosome)
        rev (seq record) - str of reverse complement sequence (should be one per chromosome)
    """

    print("reading fa", build)

    fa_dict = {
        "hs1": "/wynton/home/ahituv/fongsl/dna/hs1/hs1.fa.gz",
        "hg38": "/wynton/group/databases/goldenPath/hg38/bigZips/hg38.fa.gz",
        "saccer3": "/wynton/home/ahituv/fongsl/dna/sacCer3/sacCer3.fa.gz"
    }

    FA = fa_dict[build.lower()]

    forward, reverse = {}, {}

    # 2 open
    with gzip.open(FA, "rt") as handle:

        ### COULD PARALLEL HERE###
        
        # 3  read using simpleFastaParser
        for val in SimpleFastaParser(handle):
            start = timer()
            seq_id, seq = val

            f_dict = getKmerCount(key_pool, keysize, window_size, seq)
            forward[seq_id] = f_dict

            # 4
            seq_rev = str(Seq(seq).reverse_complement())
            r_dict = getKmerCount(key_pool, keysize, window_size, seq_rev)
            reverse[seq_id] = r_dict

            print("done w seq", seq_id, "in", timer()-start)

    return forward, reverse


def getKmerCount(key_pool, keysize, windowsize, sequence):
    """
    get 1 kmer; break up sequence into equally sized kmers w sliding window of with windowsize, stepsize

    required packages
        - numpy

    inputs
        key_pool (list) - list of keys key to filter sequences by 
        keysize (int) - size of key
        windowsize (int) - windowsize to make for sequence
        sequence (str) - sequence to break into kmers


    method
        1. establish all the positions for the sequence, 
            1.1 make a dictionary of dictionaries to count kmer keys, 
        2. iterate through start positions, get kmer, kmer key
        3. if kmer does not have an N (i.e. absent of A,C,T,G nucleotide base)
            in sequence and key in keypool, add to dictionary of dictionaries

    return 
        key_dict (dict) - dict of keys + dict of kmers and counts that match key_pool in input sequence

    """

    # 1
    n_kmers, start = (len(sequence) - windowsize + 1), 0

    #1.1 make dictionary of dictionaries w/ counters
    key_dict = {}  
    
    for k in key_pool:
        key_dict[k] = Counter()  # counter for each prefix key

    # 2
    for start in np.arange(n_kmers):

        # get kmer from sequence
        kmer = sequence[start:start+windowsize].upper()
        kmer_key = kmer[:keysize]  # get sequence prefix key

        # 3 append sequence windows within range of possible windows
        if kmer_key in key_pool:  
            
            # double check that kmer key prefix is in the key dictionary. 
            # If not, add. 
            
            # either way, invoke the key_counter variable to count that prefix. 
            if kmer_key not in key_dict:
                key_dict[kmer_key] = Counter()
                key_counter = key_dict[kmer_key]

            else:
                key_counter = key_dict[kmer_key]
            
            # count the prefix. 
            if "N" not in kmer.upper():
                key_counter[kmer] += 1

    return key_dict


def countKmersUniverse(forward_dict, reverse_dict, keysize):
    """
    add kmer counts across chromosomes per key

    input
        forward_dict (dict) - dict of forward kmers + counts in fasta sequence
        reverse_dict (dict) - dict of reverse kmers + counts in fasta sequence
        keysize (int) - size of prefix key identifier sequence

    method
        1. make key_dict to build kmer-universe dict while parsing kmers. Reduces memory burden. 
            1.1 run forward and reverse dictionaries sequentially
        2. per forward/reverse dictionary (keys = chr numbers, values = dictionary of dictionaries). Need only values here. 
            2.1 unpack dictionary of key, count_dict for each chromosome
            2.2 if key sequence is not in key dictionary, add as Counter
            2.3 iterate though key-related kmers, keep only kmer-snippets linked to key, kmer, and their counts,
                increasing the kmer count for a kmer while iterating across forward and reverse, across chromosomes.  

    return
        key_dict (dictionary) - genome universe w/ kmer counts.

    """

    key_dict = {}  # make dictionary of dictionaries w/ counters

    start = timer()

    # 2 count forward/reverse items in input dicts, 
    # where keys are chromosomes, values are 
    # dictionaries of dictionaries of keys, full kmers, and kmer counts
    
    seqs = [forward_dict, reverse_dict]

    for seq_dict in seqs:
        for chr_dict in seq_dict.values():

            # 2.1
            for key, count_dict in chr_dict.items():

                # 2.2 check that key is in dictionary
                if key not in key_dict:
                    key_dict[key] = Counter()

                # get key-specific kmer counter.
                key_counter = key_dict[key]

                # 2.3
                for seq, count in count_dict.items():
                    
                    # cut sequence into kmer snippet
                    kmer_snippet = seq[keysize:]
                    key_counter[kmer_snippet] += count
    end = timer()

    print("Counting kmer instances:", end-start)

    return dict(key_dict)


def writeDictionary(path, windowsize, result_dict, writenull):
    """
    write kmer universe dictionary of dictionaries to outfiles, write nullomers

    input
        path (str) - path to write outfile
        windowsize (int) - size of window used for kmers
        results_dict (dictionary) - dictionary of dictionaries, 
            # kmer universe w/ prefix keys, kmer sequences, 
            # and counts of kmer occurrences 
        writenull (bool) - write nullomer file or don't

    method
        1. create the outfile, sequence universe
        2. write key, value as comma-separated str to kmer outfile
            2.1. zip the file
        3. write left-over sequences to nullomer file
            3.1 zip the file

    return
        written_keys(list) - list of prefix keys
    """    
    written_keys = []
    
    for key, kmer_dict in result_dict.items():
        print("\n Writing dict", key)
        # 1
        kmer_file = os.path.join(path, f"{key}.kmers.csv")
        null_file = os.path.join(path, f"{key}.nullomers.csv")

        ### FUTURE SARAH - YOU WILL HAVE TO FIND A BETTER WAY TO MAKE PREFIX KEYS FOR LARGER SEQUENCES ###

        if writenull is True:
            # sometimes the universe is too big to make (e.g. 4**18 = 6.9e10 sequence)
            # set writenull to False
            
            # generate all potential suffix sequences 
            # of the window size minus the prefix key size.
           
            universe = makeKeys(windowsize-len(key))

        # 2
        with open(kmer_file, "w") as kmer:

            # write kmers first
            for value_seq, count in kmer_dict.items():

                # put full prefix+suffix sequence back together
                seq = key+value_seq

                kmer.write(f"{seq},{count}\n")

                if writenull is True:
                    
                    # remove suffix from potential kmer prefix universe
                    universe.remove(value_seq)

            kmer.close()
            written_keys.append(key)
            
            #2.1 zip the file
            os.system(f"gzip {kmer_file}")

        # 3
        if writenull is True:
            
            with open(null_file, "w") as nullomer:
                
                for seq in universe:
                    nullomer.write(f"{key+seq}\n")
                    
                nullomer.close()
            
            #3.1 zip the file
            os.system(f"gzip {null_file}")
            
    return written_keys


def readKeyLog(path):
    """
    read chr.log and determine whether chromosome has already been summed into file dataset

   input
        chr_num (str) - chromosome number
        path (str) - path to directory to write log

    method
        1. make the log file
            1.1 if log file does not exist, return False (none of the chromosomes have been run) 
        2. open existing log file. 
        3. append chr_num to list
        4. check whether input chr_num is in list:
            4.1 if yes, return True (chromosome has been run)
            4.2 if no, return False (chromos dome has not been run)

    return
        empty set
        chr_runlist (set) - if key has been added to final nullomer count
    """

    # 1
    out = os.path.join(path, "chr.log")

    # 2
    runlist = set()
    if os.path.exists(out) is True:

        with open(out, "r") as chrlog:
            # 3
            for line in chrlog.readlines():
                runlist.add(str(line.split("\n")[0]))

    return runlist


def writeKeyLog(key_pool, path):
    """
    write log of chrs written to the all mer

    input
        chr_num (str) - chromosome number
        path (str) - path to directory to write log

    method
        1. make the log file
        2. open the log file
        3. append the chr_num to the log file
        4. close the file
    """
    # 1
    out = os.path.join(path, "chr.log")

    # 2
    with open(out, "a") as chrlog:
        for key in key_pool:
            # 3
            chrlog.write(f"{key}\n")
    # 4
    chrlog.close()

# MAIN

def main(argv):
    """
    count number of kmer occurrences per chromosome, genome, write output to file

    method
        0. read array, get key_pool from array, make outdirectory

        1. read log, check if you've run any prefix keys already

        2. if prefix key has not been run (set), run it

        3. get genome sequence, forward and reverse as dictionaries 
            for each chromosome, each prefix key, kmers, and counts. 
            (Should be 25 items corresponding to 25 chromosomes)
            3.1 scan genome for kmer windows, keep only windows with prefix keys in the 5'
            
        4. combine prefix key, kmer counts across forward and reverse and across chromosomes. 

        5. write kmer prefix keys to log. 
    """

    # 0
    # e.g. prefix KEY_POOL = ["AAAA", "TTTT"]
    KEY_POOL, PATH, BUILD, WINDOW_SIZE, KEYSIZE, WRITENULL = array_reader(ARRAY, JOBNUM)  

    OUTDIR = makeOutDir(PATH)
    
    # 1
    run_already = readKeyLog(OUTDIR)  # check whether prefix key has been run already
    
    KEY_POOL_ = list(set(KEY_POOL).difference(run_already))  # remove run already set w set difference function, turn into list
    
    # 2
    if len(KEY_POOL_)>0:  
        print("get nullomers space from kmers", KEY_POOL_)

        # 3 get forward, reverse complement genome sequences, filtering on prefix key
        FWD, REV = extractFaSeq(PATH, KEY_POOL_, KEYSIZE, WINDOW_SIZE, BUILD)

    # 4 count knownmers in kmer universe
    key_universe = countKmersUniverse(FWD, REV, KEYSIZE)

    #5 write kmer dictionaries summed across sequences per prefix key identity. 
    written_keys = writeDictionary(OUTDIR, WINDOW_SIZE, key_universe, WRITENULL)
    
    # write the prefix key to the out dir
    writeKeyLog(written_keys, OUTDIR)
    
if __name__ == '__main__':
    main(sys.argv[1:])    