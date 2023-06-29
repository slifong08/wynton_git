import argparse
import config_readwrite as crw
import os, sys

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("build", type=str, help='genome build')
arg_parser.add_argument("length", type=str, help='length')


args = arg_parser.parse_args()

BUILD, LEN = args.build, args.length

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
    
    return str(kmer_key[int(kmer_len)])

# read config
config_tag = f"config.{BUILD}.pool.ini"
config, cfn = crw.read_config(os.path.join(os.getcwd(), config_tag))

# add section
crw.check_section(config, f"{LEN}mer")

# write
SRC_CHR = f"/wynton/home/ahituv/fongsl/dna/{BUILD}/chromosomes/"

PATH = os.path.join(os.path.dirname(os.path.dirname(SRC_CHR)), f"kmers/{LEN}mers/")
print(PATH)

RESULTS=f"/wynton/home/ahituv/fongsl/nullomers/results/{LEN}mers/"

NULLOMERS=os.path.join(PATH, f"ALL.{LEN}mers-nullomers.csv.gz")

KMERS=os.path.join(PATH, f"ALL.{LEN}mers-kmers.csv.gz")

KEY = getKeySize(LEN)
ARRAY = f"/wynton/home/ahituv/fongsl/nullomers/bin-generate/arrays/array-{BUILD}.{LEN}mers.{KEY}.tsv"

config[f"{LEN}mer"]["src_chr"] = SRC_CHR
config[f"{LEN}mer"]["path"] = PATH
config[f"{LEN}mer"]["nullomers"] = NULLOMERS
config[f"{LEN}mer"]["kmers"] = KMERS
config[f"{LEN}mer"]["results"] = RESULTS
config[f"{LEN}mer"]["keysize"] = KEY
config[f"{LEN}mer"]["kmerlen"] = LEN
config[f"{LEN}mer"]["array"] = ARRAY

crw.write_config(config, cfn)
print("made config - ", cfn)

# make paths
if os.path.exists(os.path.dirname(PATH)) is False:
    os.mkdir(os.path.dirname(PATH))

if os.path.exists(PATH) is False:
    os.mkdir(PATH)
    
if os.path.exists(RESULTS) is False:
    os.mkdir(RESULTS)
             