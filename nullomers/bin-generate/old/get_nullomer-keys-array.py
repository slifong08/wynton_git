import argparse
from collections import Counter
import glob
import gzip
from itertools import product
import numpy as np
import os, sys
import pyarrow as pa
from timeit import default_timer as timer

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("job_num", type=str, help='job number to look up in array.tsv')
arg_parser.add_argument("array", type=str, help='ARRAY file of job numbers, chromosome ids')
arg_parser.add_argument("directory", type=str, help='directory where chromosom.fa lives')
arg_parser.add_argument("length", type=int, help='kmer length')
arg_parser.add_argument("keysize", type=int, help='kmer key length for storing info')
arg_parser.add_argument("writenull", type=bool, help='write nullomers to file? (bool). Not recommended for long sequences. ')

args = arg_parser.parse_args()
JOBNUM, ARRAY, PATH, WINDOW_SIZE, KEYSIZE, WRITENULL=  args.job_num, args.array, args.directory, args.length, args.keysize, args.writenull

#JOBNUM, ARRAY, PATH, WINDOW_SIZE, KEYSIZE = 3, "/wynton/home/ahituv/fongsl/nullomers/bin-generate/arrays/array-hs1.23mer.tsv","/wynton/home/ahituv/fongsl/dna/hs1/", 23, 5

OUTDIR = os.path.join(PATH, "kmers", f"{WINDOW_SIZE}mers")

# FUNCTIONS

def array_reader(array, job_number):
    """
    read array, return values matching job number
    """
    
    with open(array, "r") as reader:
        for line in reader:
            num, keyfile, key = line.strip("\n").split("\t")
            if int(num) == int(job_number):
                KEYFILE, KEY = keyfile, key
            
    return KEYFILE, KEY


def makeChrList():
    """
    return list of chromosome numbers
    """
    chr_list = list(np.arange(1,23))
    chr_list.extend(["X", "Y"])
    
    return chr_list


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

def countKeymerSpace(windowsize, path, key):
    """
    input
        universe_kmer - dictionary with the kmer universe    
        windowsize(int) - size of the kmer universe
        path(str) - path to dataframes
        key (str) - sequence key that kmers are split on

    method
        1. concatenate kmers w/ same key across chromosomes. 

        3. per key, read sub-kmer (i.e. kmer = key + sub-kmer) counts (val) associated with that key

        4. if key is in keymer universe, than sum add the new count to the sum. 

    """
    start = timer()

    #1
    query = os.path.join(path, f"chr*.{key}.csv.gz")
    out = os.path.join(path, f"{key}.csv")
    
    os.system(f"zcat {query} > {out}")# | rm {query}")
    
    #3
    key_universe=Counter()
    
    with open(out, "r") as reader:
        for line in reader.readlines():
            
            value_seq, count = line.strip("\n").split(",")
            key_universe[value_seq] +=int(count)
            
        reader.close()
        print(f"cat {key}", out, timer()-start)

        
    return dict(key_universe), query
               
    
    
def writeDictionary(key, path, windowsize, result_dict, writenull):
    
    """
    write kmer universe dictionary to outfile, write nullomers
    
    input
        
        path (str) - path to write outfile
        windowsize (int) - size of window used for kmers
        null_file (str) - from read_array
        results_dict (dictionary) - kmer universe w/ counts of kmer occurrences 
    
    method
        1. create the outfile, sequence universe
        2. write key, value as comma-separated str to kmer outfile
        3. write left-over sequences to nullomer file
        
    return 
        out_file (str) - written file name
    
    """
    
    #1
    kmer_file = os.path.join(path, f"{key}.csv")
    null_file = os.path.join(path, f"{key}.nullomers.csv")
    
    ### FUTURE SARAH - YOU WILL HAVE TO FIND A BETTER WAY TO MAKE KEYS FOR LARGER SEQUENCES ###
    
    if writenull is True:
        universe = makeKeys(windowsize-len(key))
    
    #2
    with open(kmer_file, "w") as kmer:
        
        # write kmers first
        for value_seq in result_dict.keys():
            
            count = result_dict[value_seq]

            # full sequence back together
            seq = key+value_seq 
           
            kmer.write(f"{seq},{count}\n")

            if writenull is True:
                # remove key from kmer universe
                universe.remove(value_seq)
            
        kmer.close()
    
    #3
    if writenull is True:
        with open(null_file, "w") as nullomer:
            for seq in universe:
                nullomer.write(f"{key+seq}\n")
            nullomer.close()
    
    return kmer_file, null_file


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
        chr_runlist (set) - if key has been added to final nullomer count
    """
    
    #1
    out = os.path.join(path, "chr.log")

    #2
    runlist=set()
    if os.path.exists(out) is True:
    
        with open(out, "r") as chrlog:
        #3
            for line in chrlog.readlines():
                runlist.add(str(line.split("\n")[0]))
        
        
    return runlist
       
            
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
    out = os.path.join(path, "chr.log")
    
    #2
    with open(out, "a") as chrlog:
        
        #3
        chrlog.write(f"{chr_num}\n")
    
    #4
    chrlog.close()
    
# MAIN

def main(argv):
    """
    count number of kmer occurrences per chromosome, write output to file
    
    method
        0. read array, get keyfile, key. 
            
        1. read log, check if you've run any keys already
        
        3. if key has been run (set), take difference of key_set from run_already set. 
            Else, take difference from empty set (i.e. difference is nothing and full kmer_key set gets run) 
        4. for the kmer key, 
            4.1 Update kmer-key universe dict values w/ sum the number of instances across chromosomes
            4.2 Write the universe dictionary for every kmer key. 
            4.3 Write kmer key to log file (this line is inside the countKeymerSpace function). 
                So, only when kmer key is complete will universe be written and chrlog be updated. 
            return key universe, query for keys
        5. write key_universe dictionary to file. 
        6. write the key to the chrlog
        7. delete chr based files
    """
    

    #0
    KEYFILE, KEY = array_reader(ARRAY, JOBNUM)
    
    #1
    run_already = readChrLog(OUTDIR) # check whether key has been run already

    #2 get set of kmer keyes based on keysize 
    key_set = makeKeys(KEYSIZE)

    #3
    if KEY not in run_already:  # remove run already set. 
        
        print("get nullomers space from kmers", KEY)
        
        #4 write kmer dictionaries summed across sequences w key identity. 
        key_universe, query = countKeymerSpace(WINDOW_SIZE, OUTDIR, KEY)
        
        kmerfile, nullomerfile = writeDictionary( KEY, OUTDIR, WINDOW_SIZE, key_universe, WRITENULL)
        
        # count kmer file lines.
        linecount = sum(1 for i in open(kmerfile, "r"))
        
        # if data is written
        if linecount > 0:
            
            # write the key to the out dir
            writeChrLog(KEY, OUTDIR)

            os.system(f"rm {query}")
        else:
            print(KEY, "fail")
        

if __name__ == '__main__':
    main(sys.argv[1:])    