import argparse
import config_readwrite as crw
import glob
from itertools import product
import os, sys

"""
make an array of tags for combining nullomers, mutagenesis of nulloemrs

input
    config (str) - config
    build (str) - genome builc
    length (int) - kmer length
    keysize (int) - size of sequence key identifier

return
    array - of job numbers, combined nullomer, and key ids written to file. For running nullomers, N order sequence job arrays
"""

# args
arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("config", type=str, help='config name')
arg_parser.add_argument("build", type=str, help='genome build')
arg_parser.add_argument("length", type=int, help='kmer length')
arg_parser.add_argument("keysize", type=int, help='size of key')

args = arg_parser.parse_args()

CONFIGTAG, BUILD, KMER_LEN, KEYSIZE = args.config, args.build, args.length, args.keysize

# functions
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

def main(argv):
    """
    method
        1. make config section variable
        2. make array file name
        3. read config, get out directory
        4. make sequence keys (e.g. AAA, AAC, AAG, etc.)
        5. make nullomer file names matching sequence keys
        6. write the job number, nullomer file names, sequence keys to array
        7. write array to config
    """
    
    #1 make config section str
    SECTION=str(KMER_LEN)+"mer"
    
    #2 make array file name
    ARRAY = os.path.join(os.getcwd(), "arrays", f"array-{BUILD}.{SECTION}.tsv")

    #3 read config
    config, cfn = crw.read(CONFIGTAG)

    PATH = config[SECTION]["path"]

    #4 get keys
    KEYS = makeKeys(KEYSIZE)
    
    #5 get files matching keys
    FILES = [os.path.join(PATH, f"ALL.{key}.{KMER_LEN}mers-nullomers.csv.gz") for key in KEYS]
    print(len(FILES))

    #6 write array file
    with open(ARRAY, "w") as writer:
        for n, file in enumerate(zip(FILES, KEYS)):
            writer.write(f"{n+1}\t{file[0]}\t{file[1]}\n")

        writer.close()

    #7 write to config
    config[SECTION]["array"] = ARRAY
    crw.write(config, cfn)
    
    
if __name__ == "__main__":
    main(sys.argv[1:])