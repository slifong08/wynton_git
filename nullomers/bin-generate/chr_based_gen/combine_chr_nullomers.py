import argparse
import gzip
from itertools import product
import numpy as np
import os, sys

from timeit import default_timer as timer


arg_parser = argparse.ArgumentParser()


arg_parser.add_argument("-d","--directory", type=str, help='directory where chromosom.fa lives')
arg_parser.add_argument("-l","--length", type=int, help='kmer length')

args = arg_parser.parse_args()

PATH, WINDOW_SIZE = args.directory, args.length

OUTDIR = os.path.join(PATH, "kmers", f"{WINDOW_SIZE}mers")
if os.path.exists(OUTDIR) is False:
    os.mkdir(OUTDIR) 

# FUNCTIONS


def makeChrList():
    """
    return list of chromosome numbers
    """
    chr_list = list(np.arange(1,23))
    chr_list.extend(["X", "Y"])
    
    return chr_list

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
    
    return kmer_universe

def writeDictionary(chr_num, path, windowsize, result_dict):
    """
    write kmer universe dictionary to outfile, write nullomers
    
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
    out_file = os.path.join(path, f"{chr_num}.{windowsize}mers.csv.gz")
    null_file = os.path.join(path, f"{chr_num}.{windowsize}mers-nullomers.csv.gz")
    
    #2
    nullomers = gzip.open(null_file, "wt")
    with gzip.open(out_file, "wt") as writer:
        
        for key, value in result_dict.items():
            writer.write(f"{key},{value}\n")

            if int(value) == 0: # if nullomer
                nullomers.write(f"{key},{value}\n")
    #3
    writer.close()
    nullomers.close()
    
    print("\n\nwrote", out_file)
    
    return out_file


def checkUniverseFile(windowsize, path):
    """
    return universe file if written
    
    input
        windowsize (str) - length of kmer
        
    method
        1. check whether file exists
        2. if not - make universe of kmers
        3. if exists - read kmer universe into dictionary
        
        
    return
        universe_kmers (dict) - dictionary of kmers, summed counts across chromosomes. 
    """
    
    # 1
    outfile = os.path.join(path, f"ALL.{windowsize}mers.csv.gz")
    
    #2
    if os.path.exists(outfile) is False:
        universe_kmers = makeKmerUniverse(windowsize)  ## make 
    #3
    else:
        universe_kmers ={}
        with gzip.open(outfile, 'rt') as universe:
            for line in universe:
                kmer, count = line.strip("\n").split(",")
                universe_kmers[kmer]=int(count)
                
    return universe_kmers
        
    
def readChrLog(chr_num, path):
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
            4.2 if no, return False (chromosome has not been run)
    
    return
        False - if chromosome has not been added to final nullomer count
        True - if chromosome has been added to final nullomer count
    """
    
    #1
    out = os.path.join(path, "chr.log")

    #1.1
    if os.path.exists(out) is False:
        return False
    #2
    else:
        chr_runlist =[]
        with open(out, "r") as chrlog:
            #3
            for line in chrlog.readlines():
                chr_runlist.append(str(line.split("\n")[0]))
        
        if str(chr_num) in chr_runlist:
            return True
        
        else:
            return False
            
            
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
        1. make a list of chromosomes to query
        2. make the kmer universe dictionary
        3. per chromosome, read the chromosome log, check if you've run this chromosome already
            3.1. if chromosome has not been run (False), iterate through chromosomes, opening csv
        4. for each kmer, update kmer universe dict values w/ sum the number of instances across chromosomes
        5. write the universe dictionary after every chromosome
        6. after writing dictionary, write chromosome to log file. Thus, only when chromosome is complete will universe be written and chrlog be updated. 
        
        
    """
    
    #1
    chr_list = makeChrList()
    #2
    universe_kmers = checkUniverseFile(WINDOW_SIZE, OUTDIR)

    #3
    for CHR_NUM in chr_list:
        run_already = readChrLog(CHR_NUM, OUTDIR) # check whether chromosome has been run already
        print(CHR_NUM, run_already)

        #3.1
        if run_already == False:
            
            #4
            F = os.path.join(OUTDIR, f"{CHR_NUM}.{WINDOW_SIZE}mers.csv.gz")
            
            with  gzip.open(F,'rt') as reader:
                for line in reader:
                    key, val = line.strip("\n").split(",")
                    universe_kmers[key] += int(val)
            #5
            outcsv = writeDictionary("ALL", OUTDIR, WINDOW_SIZE, universe_kmers)
            

            #6 write chromosome number to chr.log  
            writeChrLog(CHR_NUM, OUTDIR)
            
            os.chdir(OUTDIR)
            
            # concat KMERS
            os.system(f"cat ALL.*.{}mers.csv > ALL.{WINDOW_SIZE}.csv && gzip ALL.{WINDOW_SIZE}.csv")
            
            # concat Nullomers
            os.system(f"cat ALL.*.{}mers-nullomers.csv > ALL.{WINDOW_SIZE}-nullomers.csv && gzip ALL.{WINDOW_SIZE}-nullomers.csv"
            
        else:
            continue

if __name__ == '__main__':
    main(sys.argv[1:])