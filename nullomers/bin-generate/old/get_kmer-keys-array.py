#!/usr/bin/env python
# coding: utf-8
"""
sarahfong

Summary
    - count kmers per chromosome
    - split kmer into multiple files where the first N letters of the kmer is the "key"
    - one file per key will be written w/ counts of sequences associated w/ that key + value

Example
    Say we find an 11mer, GCGTACGTACG, that appears 15025 times on chrN. The key length is 4. 
        The key is the first 4 letters "GCGT". 
        The value, or other 7 letters, "ACGTACG" will be written to chrN.GCGT.csv
        The 11mer counts for GCGTACGTACG, will be written to chrN.GCGT.csv

        The resulting file will look like this:

        ./chrN.GCGT.csv
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

arg_parser.add_argument("job_num", type=str, help='job number to look up in array.tsv')
arg_parser.add_argument("array", type=str, help='ARRAY file of job numbers, key ids')
arg_parser.add_argument("directory", type=str, help='directory where genome.fa lives')
arg_parser.add_argument("length", type=int, help='kmer length')
arg_parser.add_argument("keysize", type=int, help='kmer key length for storing info')

args = arg_parser.parse_args()

JOB_NUM, ARRAY, PATH, WINDOW_SIZE, KEYSIZE = args.job_num, args.array, args.directory, args.length, args.keysize

"""
JOB_NUM=13 # 25 = chrY, 13=chr22
ARRAY = "/wynton/home/ahituv/fongsl/nullomers/bin-generate/arrays/chr_fa_array-hs1.tsv"
PATH = "/wynton/home/ahituv/fongsl/dna/hs1/"
WINDOW_SIZE = 12
KEYSIZE = 5
"""

# make outdirectory for
OUTDIR = os.path.join(PATH, "kmers", f"{WINDOW_SIZE}mers")

print("OUTDIR", OUTDIR) 
if os.path.exists(OUTDIR) is False:
    if os.path.exists(os.path.dirname(OUTDIR)) is False:
        
        os.mkdir(os.path.dirname(OUTDIR)
    os.mkdir(OUTDIR) 


###
# functions
###

def array_reader(array, job_number):
    """
    read array, return value matching job number
    """
    
    with open(array, "r") as reader:
        for line in reader:
            num, key = line.strip("\n").split("\t")
            if int(num) == int(job_number):
                CHR_NUM = key
            
    return CHR_NUM

# ## extract fa sequence
def extractFaSeq(chr_num, path):
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

    print("reading fa", str(chr_num) +".fa.gz")

    #1
    WYNTON_PATH = "/wynton/group/databases/goldenPath/hg38/chromosomes"
    fa = os.path.join(WYNTON_PATH, f"{chr_num}.fa.gz")

    #2 open
    with gzip.open(fa, "rt") as handle:

        #3  read using simpleFastaParser
        for val in SimpleFastaParser(handle):
            seq_id, seq = val
            print(seq_id)

            #4
            rev = str(Seq(seq).reverse_complement())

            print("sequence size forward, reverse", len(seq), len(rev))
    return seq, rev

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


def getOneKmer(sequence, windowsize, start):

    """
    get 1 kmer; break up sequence into equally sized kmers w sliding window of with windowsize, stepsize
    
    required packages
        - numpy
        
    inputs
        sequence (str) - sequence to break into kmers
        windowsize (int) - windowsize to make for sequence
        start (int) - index to get sequence

        
    method
        1. get kmer
        2. if kmer does not have N in sequence
    
    return 
        kmer (str) - windowsize from input sequence
        
    """
    #1
    kmer = sequence[start:start+windowsize].upper()

    #2
    if "N" not in kmer.upper(): # append sequence windows within range of possible windows

        return kmer
    
    else:
        return None
    
# ## count kmers in universe

def countKmersUniverse(kmer_list, keysize):
    
    """
    add kmer counts to chromosome universe dictionary:

    input
        kmer_list (list) - list of kmers in fasta sequence

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

    universe_dict = {}


    start = timer()
    
    
    #2 count items in dictionary
    for i in kmer_list:

        if i is not None:

            key, value = i[:keysize], i[keysize:] # split string into key, value

            if key not in universe_dict: # if key is not in universe dictionary
                universe_dict[key]={}

            #2.1
            keydict = universe_dict[key]

            #2.2 handle cases when value is not in dictionary
            if value not in keydict:
                keydict[value] = 1
                
            else:
                keydict[value]+=1
               
        else:
            continue
            
    end = timer()

    print("Counting kmer instances:", end-start, "\n Writing dict")
    
    return universe_dict


# ## write dictionary

def writeDictionary(chr_num, outpath, windowsize, result_dict, keyset):
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
    #print("writing keys", keyset)
    
    #1
    for key, value in result_dict.items():
        log_key = f"{chr_num}.{key}"
        
        # make the outfile for each key 
        out_file = os.path.join(outpath, f"{log_key}.csv")
        
        if log_key in keyset:

            # write the values
            with open(out_file, "w") as writer:
                for secondkey, count in value.items():

                    writer.write(f"{secondkey},{count}\n")
            #3
            writer.close()

            #4
            os.system(f"gzip {out_file}")

            #5 write outfile to log
            writeChrLog(log_key, outpath)
        
    query = os.path.join(outpath, f"{chr_num}.*.csv")
    
    return query


def readChrLog(path):
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
    
    #1
    out = os.path.join(path, "kmer.log")

    #2
    runset=set()
    if os.path.exists(out) is True:
    
        with open(out, "r") as chrlog:
        #3
            for line in chrlog.readlines():
                runset.add(str(line.split("\n")[0]))
        
    print("N chromosomes have been run already", len(runset))
    return runset
       
            
def writeChrLog(chr_num, path):
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
    #1
    out = os.path.join(path, "kmer.log")
    
    #2
    with open(out, "a") as chrlog:
        
        #3
        chrlog.write(f"{chr_num}\n")
    
    #4
    chrlog.close()
    
    #print(f"wrote {chr_num} to kmer.log", out)
    

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
        3. start the timer. 
        4. set up to retrieve partial commands
        5. pool forward, reverse kmers from sequence
        6. extend and count kmers in genome
        7. write dictionary to .csv

    return
        outcsv (str) - csvs of kmer universe, separated by keys
    
    """
    
    #0
    CHR_NUM = array_reader(ARRAY, JOB_NUM)
    
    #1 get chr sequence, reverse complement 
    SEQ, REV = extractFaSeq(CHR_NUM, PATH)


    #2 get kmer key set, universe dict
    KEY_SET = makeKeys(KEYSIZE)
    
    CHR_KEYS = set([f"{CHR_NUM}.{key}" for key in KEY_SET])
    
    # 
    run_already = readChrLog(OUTDIR) # check whether key has been run already
    

    WRITE_KEYS = set(CHR_KEYS).difference(run_already) 


    #3 calculate the number of starting positions 
    STARTS = list(range(len(SEQ) - WINDOW_SIZE + 1))

    #4 setup partial command to retrieve kmers
    partial_seq = partial(getOneKmer, SEQ, WINDOW_SIZE)
    partial_rev = partial(getOneKmer, REV, WINDOW_SIZE)

    #5 get all kmers from chromosome
    result_multi = [partial_seq(start) for start in STARTS]
    result_multi_r = [partial_rev(start) for start in STARTS]

    #6 count knownmers in kmer universe
    result_multi.extend(result_multi_r) # make 2 lists into 1 big list

    # count items, add to dictionary
    genome_kmers = countKmersUniverse(result_multi, KEYSIZE)

    #7 write dictionary
    outcsv = writeDictionary(CHR_NUM, OUTDIR, WINDOW_SIZE, genome_kmers, WRITE_KEYS)

    print(outcsv)
    
if __name__ == '__main__':
    main(sys.argv[1:])