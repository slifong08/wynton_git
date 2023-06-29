#!/usr/bin/env python
# coding: utf-8

# In[4]:


import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from itertools import product
import numpy as np
import os, sys

from multiprocessing import Pool, cpu_count

from timeit import default_timer as timer
from functools import partial


# # parse arguments

# In[5]:


arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("-c","--chromosome", type=str, help='chromosome number')
arg_parser.add_argument("-d","--directory", type=str, help='directory where chromosom.fa lives')
arg_parser.add_argument("-l","--length", type=int, help='kmer length')

args = arg_parser.parse_args()

CHR_NUM, PATH, WINDOW_SIZE = args.chromosome, args.directory, args.length

# make outdirectory for
OUTDIR = os.path.join(PATH, "kmers", f"{WINDOW_SIZE}mers")
if os.path.exists(OUTDIR) is False:
    os.mkdir(OUTDIR) 

# parse if chromosome number should be x, y from array input. 
if CHR_NUM == "23":
    CHR_NUM = "X"
elif CHR_NUM == "24":
    CHR_NUM = "Y"


# # functions
# ## extract fa sequence

# In[ ]:


def extractFaSeq(chr_num, path):
    """
    get chr-specific sequence, reverse complement
    
    require
        Bio.SeqIO.FastaIO - SimpleFastaParser
    
    input
        chr_num (str) - chromosome number e.g. 1,2,3,X,Y
        path (str) - path to .fa files
        
    method
        1. get fa file, check if fa file is zipped
        2. open fa file
        3. get sequence w/ simpleFastaParser
        4. get reverse complement by turning seq into Seq object
        5. rezip .fa
    
    return 
        seq (seq record) - str of sequence (should be one per chromosome)
        rev (seq record) - str of reverse complement sequence (should be one per chromosome)
    """
    start = timer()
    print("extracting", chr_num+".fa")
    
    #1
    fa = os.path.join(path, f"chr{chr_num}.fa")
    zipped = fa+".gz"
    
    if os.path.exists(zipped) is True:
        print('unzipping')
        os.system(f'gunzip {zipped}')

    #2 open
    with open(fa) as handle:

        #3  read using simpleFastaParser
        for val in SimpleFastaParser(handle):
            seq_id, seq = val

            #4
            rev = str(Seq(seq).reverse_complement())
    #5 rezip
    if os.path.exists(fa):
        print("rezipping")
        os.system(f'gzip {fa}')
    
    print("done extracting", timer() - start)
    
    return seq, rev


# ## make kmer universe

# In[ ]:


def makeKmerUniverse(LEN):
    """
    input
        len (int) - length of kmer to create
    
    require
        itertools.product
        python list comprehension, join

    method
        itertools.product kmers -> join kmers into str -> make list

    return
        kmer_universe (dictionary) - list of strings. 
    """
    kmer_universe = {}
    v=0
    for item in product("ACGT", repeat=LEN):
        seq="".join(item)
        kmer_universe[seq]=0
        #v+=1
    
    return kmer_universe


# ## retreive one kmer
# 
# ### !Major question! What to do about N's

# In[ ]:


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


# ## pool kmer retreival

# In[ ]:


def poolKmer(partial_cmd, starts):
    """
    time pooled retreival of kmers from mapping starts. 
    
    require
        multiprocessing Pool
    
    input
        partial_cmd (partial object) - command to run in parallel
        starts (list) - list of starting positions to retreive kmers from in sequences
    
    method 
        1. start timer
        2. pool kmer retreival
        3. end timer
        
    return
        result_multi (list) - kmers from the parallel process
    """
    start = timer()
    n = cpu_count()
    with Pool(n) as p:

        result_multi = p.map(partial_cmd, starts) # process forward

    end = timer()
    print("Multiprocessing kmers: ", end - start)
    
    return result_multi


# ## count kmers in universe

# In[ ]:


def countKmersUniverse(kmer_list, universe_dict):
    
    """
    add kmer counts to universe dictionary
    
    input
        kmer_list (list) - list of kmers in fasta sequence
        universe_dict (dictionary) - dictionary of kmer space. Keys are kmers, values are counts
    method
        1. add kmer counts to dictionary, removing NoneType instances. 
        2. print timer
    
    return
        universe_dict (dictionary) - UPDATED w/ kmer counts. 
        
    """
    # count knownmers in kmer universe
    start = timer()

    # count items in dictionary
    for i in kmer_list:
        if i is not None:
            universe_dict[i]+=1
        else:
            continue
            
    end = timer()

    print("Counting kmer instances:", end-start, "\n Writing dict")
    
    return universe_dict


# ## write dictionary

# In[ ]:


def writeDictionary(chr_num, path, windowsize, result_dict):
    """
    write kmer universe dictionary to outfile
    
    input
        chr_num (str) - chr number
        path (str) - path to write outfile
        windowsize (int) - size of window used for kmers
        results_dict (dictionary) - kmer universe w/ counts of kmer occurrences 
    
    method
        1. create the outfile
        2. write key, value as comma-separated str to outfile
        3. close the outfile
    
    return 
        out_file (str) - written file name
    
    """
    
    #1
    out_file = os.path.join(path, f"{chr_num}.{windowsize}mers.csv")
    
    #2
    with open(out_file, "w") as writer:
        for key, value in result_dict.items():
            writer.write(f"{key},{value}\n")
    #3
    writer.close()
    
    print("\n\nwrote", out_file)
    
    return out_file


# # main


# script taken from: https://stackoverflow.com/questions/63096168/how-to-apply-multiprocessing-to-a-sliding-window 
# 
#     1e6 sequences
#     
#     Single core:  0.4470503553748131
#     
#     Multiprocessing:  0.5516961812973022

# In[ ]:


def main(argv):
    
    """
    per chromosome, write kmer count frequency dictionary as csv
    
    input
        chr_num (str) - chromozome number to process kmers for
        path (str) - path to fa files, and to write
        window_size (int) - size of window to consider for kmer
    
    method
        1. get sequence forward and reverse from fa
        2. generate kmer universe
        3. make list that ranges number of possible kmers
        4. set up to retrieve partial commands
        5. pool forward, reverse kmers from sequence
        6. extend and count kmers in universe
        7. write dictionary to .csv

    return
        outcsv (str) - csv of kmer universe
    
    """
    outputfile = os.path.join(PATH, f"{CHR_NUM}.{WINDOW_SIZE}mers.csv")
    
    if os.path.exists(outputfile) is False and os.path.exists(outputfile+".gz") is False:

        #1 get chr sequence forward, reverse complement 
        SEQ, REV = extractFaSeq(CHR_NUM, PATH)

        #2 get kmer universe dict
        universe_kmers = makeKmerUniverse(WINDOW_SIZE)

        #3 calculate the number of starting positions 
        STARTS = list(range(len(SEQ) - WINDOW_SIZE - 1))

        #4 setup partial command to retrieve kmers
        partial_seq = partial(getOneKmer, SEQ, WINDOW_SIZE)
        partial_rev = partial(getOneKmer, REV, WINDOW_SIZE)

        #5 get all kmers from chromosome
        result_multi = poolKmer(partial_seq, STARTS)  # forward strand kmers 
        result_multi_r = poolKmer(partial_rev, STARTS)  # reverse strand kmers

        #6 count knownmers in kmer universe
        result_multi.extend(result_multi_r) # make 2 lists into 1 big list

        # count items, add to dictionary
        universe_kmers = countKmersUniverse(result_multi, universe_kmers)

        #7 write dictionary
        outcsv = writeDictionary(CHR_NUM, OUTDIR, WINDOW_SIZE, universe_kmers)
    else:
        print("already looked at these kmers?", outputfile) 
    

if __name__ == '__main__':
    main(sys.argv[1:])


# # need function to write and consolidate dictionaries. 
