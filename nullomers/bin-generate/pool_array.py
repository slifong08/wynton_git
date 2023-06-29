import argparse
import config_readwrite as crw
import glob
from itertools import product
import os, sys

"""
make an array of tags for combining nullomers, mutagenesis of nulloemrs

input

    path(str) - path to output dir
    build (str) - genome builc
    length (int) - kmer length
    keysize (int) - size of sequence key identifier


return
    array - of job numbers, combined nullomer, and key ids written to file. For running nullomers, N order sequence job arrays
"""


# args
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("config", type=str, help='config name')
arg_parser.add_argument("path", type=str, help='outdir')
arg_parser.add_argument("build", type=str, help='genome build')
arg_parser.add_argument("length", type=int, help='kmer length')
arg_parser.add_argument("nmuts", type=int, help='n order mutants')


args = arg_parser.parse_args()

CONFIG, PATH, BUILD, KMER_LEN, NMUTS = args.config, args.path, args.build, args.length, args.nmuts


"""
BUILD = "sacCer3"  # "rhemac10"# "hs1" #"hg38"


PATH = f"/wynton/home/ahituv/fongsl/dna/{BUILD}/"
KMER_LEN = "11"
KEYSIZE = "4"

WRITENULL = True
PATH, BUILD, KMER_LEN, KEYSIZE
"""

# functions
def getKeySize(kmer_len):
    
    """
    determine the key sequence size to split up sequences. 
    """
    kmer_key = {11:4,
                12:4,
                13:4, 
                14:4, 
                15:4,
                16:5, 
                17:6,
                18:7,
                19:7,
                20:8,
                21:8,
                22:9,
                23:9           
    }
    
    return kmer_key[kmer_len]

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

def makePools(keys, kmer_len):
    """
    make list of keys in units of 10 
    """
    
    sizenpool = {11:10,  # dictionary of kmer size and list size to write given kmer size
                 12:10,
                 13:10,
                 14:10,
                 15:10, 
                 16:10,
                 17:10,
                 18:10,
                 19:10,
                 20:10,
                 21:10, 
                 22:10, 
                 23:10      
    }
    
    key_pools, current_pool = [], []  # empty lists to collect key pools, current pool
    
    last_key = len(keys)-1  # get index of last key
    
    listsize = sizenpool[kmer_len]  # get desired list size based on kmer size
    
    
    for n, k in enumerate(keys):
        current_pool.append(k)  # append the key to the list
        
        if n>0 and n%listsize==0:  ## for every unit of listsize (except the first unit)
            key_pools.append(current_pool) # append current list to key list, 
            current_pool=[]  # reset current_list

        elif n == last_key:  # if this is the last item in the key list
            key_pools.append(current_pool)  ## append current_list to key_pool list
    
    return key_pools  # list of lists w equal number of keys


def main(argv):
    """
    method
        0. open the config, make a section for these kmers
        1. make array file name
        2. make sequence keys (e.g. AAA, AAC, AAG, etc.)
        3. make pools of keys, with 10 keys to a pool. 
        4. set WRITENULL variable
        5. create OUTDIR
        6. write the job number, nullomer file names, sequence keys to array
        7. write array to config

    """
    
    #0
    config, cfn = crw.read(CONFIG) 
    SECTION = f"{str(KMER_LEN)}mers"
    crw.check(config, SECTION)
    
    KEYSIZE= getKeySize(KMER_LEN)
    
    #1 make array file name
    ARRAY = os.path.join(os.getcwd(), "arrays", f"array-{BUILD}.{SECTION}.{KEYSIZE}mer.{NMUTS}.tsv")
    
    #2 get keys
    KEYS = makeKeys(int(KEYSIZE))
    
    #3 make pools
    KEY_POOLS = makePools(KEYS, int(KMER_LEN))
    
    #4 do not write 18mer nullomers
    WRITENULL = True if int(KMER_LEN)<18 else False
    
    
    #5 outpath
    OUTPATH = os.path.join(PATH, "kmers", f"{KMER_LEN}mers")


    #6 write array file
    with open(ARRAY, "w") as writer:
        for n, key_pool in enumerate(KEY_POOLS):
            
            pool_str = ",".join(key_pool)  # make list into str
            
            writer.write(f"{n+1}\t{pool_str}\t{OUTPATH}\t{BUILD}\t{KMER_LEN}\t{KEYSIZE}\t{WRITENULL}\t{NMUTS}\n")

        writer.close()

    #7
    config[SECTION]["array"] = ARRAY
    crw.write(config, cfn)
    
    print("wrote pool array", ARRAY)
    
if __name__ == "__main__":
    main(sys.argv[1:])

