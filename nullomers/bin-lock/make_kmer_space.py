#!/usr/bin/env python
# coding: utf-8
"""
sarahfong

Summary
    - count kmers per fasta
    - split kmer into multiple files where the first N letters of the kmer is the "key"
    - one file per key will be written w/ counts of sequences associated w/ that key + value

Example
    Say we find an 11mer, GCGTACGTACG, that appears 15025 times on chrN. The key length is 4. 
        The key is the first 4 letters "GCGT". 
        The value, or other 7 letters, "ACGTACG" will be written to GCGT.csv
        The 11mer counts for GCGTACGTACG, will be written to GCGT.csv

        The resulting file will look like this:

        ./GCGT.csv
            ACGTACG, 15025
Limitations
    - Skips sequences with "Ns"
            
"""

import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
import glob
import gzip
from itertools import product
import numpy as np
import os, sys

from timeit import default_timer as timer
from functools import partial

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("fasta", type=str, help='.fa file to make kmer space from (full path)')
arg_parser.add_argument("length", type=int, help='kmer length')
arg_parser.add_argument("keysize", type=int, help='kmer key length for storing info')


args = arg_parser.parse_args()

FASTA, WINDOW_SIZE, KEYSIZE = args.fasta,args.length, args.keysize

"""
FASTA="/wynton/home/ahituv/fongsl/nullomers/data/lock/hepg2/filtered.hepg2.cCRE.liftOver.to.hs1.fa" 
WINDOW_SIZE=23
KEYSIZE=5
"""

# get path, name from fasta file. 
PATH, NAME = os.path.split(FASTA)[0], os.path.split(FASTA)[1].strip('.fa')

# make outdirectory for results
OUTDIR = os.path.join(PATH, "kmers", f"{WINDOW_SIZE}mers")

print("OUTDIR", OUTDIR) 

# check to make sure directories exist, or else make them
if os.path.exists(OUTDIR) is False:
    if os.path.exists(os.path.dirname(OUTDIR)) is False: # make ./kmer dir first
        os.mkdir(os.path.dirname(OUTDIR))
    os.mkdir(OUTDIR) # make full dir ./kmer/{kmer}mers

###
# functions
###

def array_reader(array, job_number):
    """
    read array, return value matching job number
    """
    
    with open(array, "r") as reader:
        for line in reader:
            num, chr_num = line.strip("\n").split("\t")
            if int(num) == int(job_number):
                CHR_NUM = chr_num
            
    return CHR_NUM

# ## extract fa sequence
def extractFaSeq(fasta):
    """
    get fasta sequence, reverse complement

    require
        Bio.SeqIO.FastaIO - SimpleFastaParser
        gzip (if zipped)

    input
        chr_num (str) - chromosome number e.g. 1,2,3,X,Y
        path (str) - path to .fa files

    method
        1. make empty forward/reverse dictionaries to collect sequences 
        2. open fa file
        3. get sequence w/ simpleFastaParser
        4. get reverse complement by turning seq into Seq object

    return 
        seq (seq record) - str of sequence (should be one per chromosome)
        rev (seq record) - str of reverse complement sequence (should be one per chromosome)
    """

    print("reading fa", fasta)

    #1
    forward, reverse = {}, {}
    
    #2
    if ".gz" in fasta:
        handle = gzip.open(fasta, "rt")
    else:
        handle = open(fasta, "r")

    #3 read using simpleFastaParser 
        
    for val in SimpleFastaParser(handle):
        seq_id, seq = val
        forward[seq_id + "_+"] = seq
        #print(seq_id)

        #4 reverse complement
        rev = str(Seq(seq).reverse_complement())
        reverse[seq_id + "_-"] = seq
    
    print("sequence size forward, reverse", len(forward.keys()), len(reverse.keys()))
            
    return forward, reverse

# make key-mer universe

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

# ## retreive one kmer
# 
# ### !Major question! What to do about N's


def getOneSeqKmers(windowsize, sequence):

    """
    get 1 sequence; break up sequence into equally sized kmers w sliding window of with windowsize, stepsize of 1
    
    required packages
        - numpy
        
    inputs
        
        windowsize (int) - windowsize to make for sequence
        sequence (str) - sequence to break into kmers

        
    method
        1. quantify the number of starts per sequence, given length. 
        2. make empty kmer list, then iterate through start positions, getting kmers. 
        3. Append to kmer list if kmer does not have N in sequence
    
    return 
        kmer_list (str) - windowsize from input sequence
        
    """
    #1
    STARTS = list(range(len(sequence) - windowsize + 1))
    
    #2
    kmer_list = []
    
    for start in STARTS:
        kmer = sequence[start:start+windowsize].upper()

        #3
        if "N" not in kmer.upper(): # append sequence windows within range of possible windows
            kmer_list.append(kmer)
    
        else:
            continue

    return kmer_list
    

# ## count kmers in universe

def countKmersUniverse(kmer_list, keysize):
    """
    add kmer counts to chromosome universe dictionary:

    input
        kmer_list (list) - list of lists of kmers in fasta sequences (forward and reverse)

    method
        1. build kmer-universe dict while parsing kmers
        2. add kmer counts to dictionary, removing NoneType instances. 
            2.1 get key-specific dictionary from universe dictionary
            2.2 if value sequence is not in key dictionary, add. 
            2.3 if value sequence has been seen, increase dictionary value count. 
        3. print timer

    return
        universe_dict (dictionary) - chromosome universe w/ kmer counts. 

    """

    # 1
    universe_dict = {}

    start = timer()

    # 2 count kmers from forward, reverse lists (kmer_list is list of lists) dictionary
    for l in kmer_list:
        for i in l:

            if i is not None:

                # split string into key, value
                key, value = i[:keysize], i[keysize:]

                # if key is not in universe dictionary, add the key
                if key not in universe_dict:
                    universe_dict[key] = {}

                # 2.1 get the key dictionary from universe dictionary
                keydict = universe_dict[key]

                # 2.2 start count when value sequence is not in dictionary
                if value not in keydict:
                    keydict[value] = 1

                # 2.3 add count of value sequences in key dictionary in universe dictionary
                else:
                    keydict[value] += 1

            else:
                continue

    end = timer()

    print("Counting kmer instances:", end-start)

    return universe_dict


# ## write dictionary

def writeDictionary(name, outpath, windowsize, result_dict): # keyset
    """
    write kmer universe dictionary to outfile

    input
        chr_num (str) - chr number or all
        path (str) - path to write outfile
        windowsize (int) - size of window used for kmers
        results_dict (dictionary) - kmer universe w/ counts of kmer occurrences 

    method
        1. create the outfile as key of keys
        2. if not already written (log_key not in keyset), write key, value as comma-separated str to outfile
        3. close the outfile
        4. compress outfile
        5. write outfile log_key to the log

    return 
        query (str) - written file name w asterisks for the key

    """
    print("Writing dict")

    # 1
    for key, value in result_dict.items():
        log_key = f"{name}.{key}"

        # 2 make the outfile for each key
        out_file = os.path.join(outpath, f"{key}.csv")

        # if key in keyset:

        # write the values
        with open(out_file, "w") as writer:
            for secondkey, count in value.items():

                writer.write(f"{secondkey},{count}\n")
        # 3
        writer.close()

        # 4
        os.system(f"gzip {out_file}")

        # 5 write outfile to log
        writeRunLog(log_key, outpath)

    query = os.path.join(outpath, f"*.csv")

    return query


def readRunLog(path):
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
        runset (set) - if key has been added to final nullomer count
    """

    # 1
    out = os.path.join(path, "kmer.log")

    # 2
    runset = set()
    if os.path.exists(out) is True:

        with open(out, "r") as runlog:
            # 3
            for line in runlog.readlines():
                runset.add(str(line.split("\n")[0]))

    print("N keys have been run already", len(runset))
    return runset


def writeRunLog(key, path):
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
    out = os.path.join(path, "kmer.log")

    # 2
    with open(out, "a") as runlog:

        # 3
        runlog.write(f"{key}\n")

    # 4
    runlog.close()

    # print(f"wrote {chr_num} to kmer.log", out)

###
# main
###

def main(argv):
   

 
    """
    per chromosome, write kmer count frequency dictionary as csv
    
    input
        chr_num (str) - chromozome number to process kmers for
        path (str) - path to fa files, and to write
        window_size (int) - size of window to consider for kmer
    
    method
        0. get chr num from job+num and array
        1. get sequence +/- from fa
        2. generate set of keys to write kmer-space to. Keys will be used to quantify the known-kmer space.
            2.1 check to make sure you haven't run this already. 
        
        3. set up to retrieve partial commands
        4. get all kmers from sequence in forward and reverse, extend to make one list of lists
        
        5. count kmers found in fasta
        7. write dictionary to .csv

    return
        outcsv (str) - csvs of kmer universe, separated by keys
    
    """

    # 0
    # NAME = array_reader(ARRAY, JOB_NUM)

    # 1 get chr sequence, reverse complement
    SEQ, REV = extractFaSeq(FASTA)

    # 2 get kmer key set, universe dict
    KEY_SET = makeKeys(KEYSIZE)

    # 2.1 
    run_already = readRunLog(OUTDIR)  # check whether key has been run already

    WRITE_KEYS = set(KEY_SET).difference(run_already)


    #3 setup partial command to retrieve kmers
    partial_seq = partial(getOneSeqKmers, WINDOW_SIZE)
    partial_rev = partial(getOneSeqKmers, WINDOW_SIZE)

    #4 get all kmers from sequences
    result_multi = [partial_seq(sequence) for sequence in SEQ.values()]
    result_multi_r = [partial_rev(sequence) for sequence in REV.values()]

    # make 2 lists into 1 big list
    result_multi.extend(result_multi_r) 

    # 5 count knownmers in kmer universe
    genome_kmers = countKmersUniverse(result_multi, KEYSIZE)

    #6 write dictionary
    outcsv = writeDictionary(NAME, OUTDIR, WINDOW_SIZE, genome_kmers) #, WRITE_KEYS)

    print(outcsv)
    
if __name__ == '__main__':
    main(sys.argv[1:])