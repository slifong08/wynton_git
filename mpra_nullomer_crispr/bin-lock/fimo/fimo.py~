import argparse

import glob
import os, sys
import pandas as pd

sys.path.append("/wynton/group/ahituv/fongsl/tools/py_")
import config_readwrite as crw

# parse args
parser = argparse.ArgumentParser()
parser.add_argument("fasta", type=str, help = "fastafile")
parser.add_argument("config", type=str, help = "config file")

args = parser.parse_args()
FA, CONFIG = args.fasta, args.config

config, cfn = crw.read(CONFIG) 

# str splitting
OUTDIR, SAMPLE_NAME = os.path.split(FA)
SAMPLE_ID= os.path.splitext(SAMPLE_NAME)[0]

# FILE CONSTANTS
FIMO_SRC = "/wynton/home/ahituv/fongsl/bin/meme-5.5.1/src/fimo"
JASPAR = "/wynton/group/ahituv/tfbs_motif/jaspar/JASPAR2022_CORE_non-redundant_pfms_meme.txt"

# write to config
section = f"FIMO"
crw.check(config, section)

config[section]["src"] = FIMO_SRC

# JASPAR file 
config[section]["jaspar"] = JASPAR

# directory to FIMO results
FIMO_RESULT_DIR = os.path.join(OUTDIR, f"fimo")
config[section]["path"] = FIMO_RESULT_DIR

crw.write(config, cfn)

# run fimo
if os.path.exists(FIMO_RESULT_DIR) is False:
    os.mkdir(FIMO_RESULT_DIR)
    
os.chdir(OUTDIR)

cmd = [
    FIMO_SRC,
    JASPAR,
    FA, 
    "--no-qvalue --text",  # do not compute FDR, do not output other datafile types 
    "&& mv" ,
    os.path.join(OUTDIR, "fimo_out"),
    os.path.join(OUTDIR, f"fimo.{SAMPLE_ID}")
    
]
print(" ".join(cmd))
os.system(" ".join(cmd))