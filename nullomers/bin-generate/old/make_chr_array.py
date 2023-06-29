"""
write array of chromosome .fa.gz files to run w/ SGE job array. 

"""
import argparse
import config_readwrite as crw
import glob
import os, sys


# args
arg_parser = argparse.ArgumentParser()

arg_parser.add_argument("path", type=str, help='path to fa.gz')
arg_parser.add_argument("build", type=str, help='genome build')
arg_parser.add_argument("config", type=str, help='config')

args = arg_parser.parse_args()

PATH, BUILD, CONFIG = args.path, args.build, args.config


# function
def getFaPath(build):
    build_dict = {
        "hg38":"/wynton/group/databases/goldenPath/hg38/chromosomes", 
        "hs1":"/wynton/home/ahituv/fongsl/dna/hs1/chromosomes", 
        "rhemac10":"/wynton/home/ahituv/fongsl/dna/rhemac10/chromosomes",
        "mm39":"/wynton/home/ahituv/fongsl/dna/mm39/chromosomes"
    }
     
    return build_dict[build]


# read config
config, cfn = crw.read(CONFIG)


# make chr array file
ARRAY = os.path.join(os.getcwd(), "arrays", f"chr_fa_array-{BUILD}.tsv")


# get FA path 
FA_PATH = getFaPath(BUILD)


# get chr FA files
FS = glob.glob(os.path.join(FA_PATH, "*.gz"))


# write array - jobnumber and chromosome.fa 
with open(ARRAY, "w") as writer:
    for i, F in enumerate(FS):
        chr_num = (os.path.split(F)[1]).strip(".fa.gz")
        writer.write(f"{i+1}\t{chr_num}\n")
    writer.close()

    
# write to cofig
SET = "ARRAY"
crw.check(config, SET)


config[SET]["FA_CHR"] = ARRAY
crw.write(config, cfn) 