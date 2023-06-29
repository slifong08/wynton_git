import argparse
import config_readwrite as crw
import os, sys

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("build", type=str, help='genome build')
arg_parser.add_argument("length", type=str, help='length')
arg_parser.add_argument("key", type=str, help='keysize')

args = arg_parser.parse_args()

BUILD, LEN, KEY = args.build, args.length, args.key

# read config
config_tag = f"config.{BUILD}.ini"
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

config[f"{LEN}mer"]["src_chr"] = SRC_CHR
config[f"{LEN}mer"]["path"] = PATH
config[f"{LEN}mer"]["nullomers"] = NULLOMERS
config[f"{LEN}mer"]["kmers"] = KMERS
config[f"{LEN}mer"]["results"] = RESULTS
config[f"{LEN}mer"]["keysize"] = KEY
config[f"{LEN}mer"]["kmerlen"] = LEN


crw.write_config(config, cfn)

# make paths
if os.path.exists(os.path.dirname(PATH)) is False:
    os.mkdir(os.path.dirname(PATH))

if os.path.exists(PATH) is False:
    os.mkdir(PATH)
    
if os.path.exists(RESULTS) is False:
    os.mkdir(RESULTS)
             