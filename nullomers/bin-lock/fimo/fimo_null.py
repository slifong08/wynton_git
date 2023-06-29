import argparse

import glob
import os, sys
import pandas as pd

sys.path.append("/wynton/home/ahituv/fongsl/tools/py_")
import config_readwrite as crw

# parse args
parser = argparse.ArgumentParser()
parser.add_argument("cellline", type=str, help = "cellline")
parser.add_argument("kmer_length", type=int, help = "length of kmer")
parser.add_argument("n_mut", type=int, help = "number of mutations to make")
parser.add_argument("config", type=str, help = "path to config")

args = parser.parse_args()
CL, KMER_LEN, NMUTS, CONFIG = args.cellline, args.kmer_length, args.n_mut,args.config

config, cfn = crw.read(CONFIG) 
"""
# dev
CL, KMER_LEN, NMUTS = "common", 14, 2
config, cfn = crw.read(os.path.join(os.getcwd(), "config.ini"))
"""

mutsection = f"{CL}.{KMER_LEN}mer.{NMUTS}mut"

# read config
OUTDIR = config[CL]["path"] 
FA = config[mutsection]["fa_fo-True"] 

# write to config
section = f"{CL}.FIMO"
crw.check(config, section)

JASPAR="/wynton/group/ahituv/tfbs_motif/jaspar/JASPAR2022_CORE_non-redundant_pfms_meme.txt"
config[section]["jaspar"] = JASPAR

FIMO_RESULT_DIR = os.path.join(OUTDIR, f"fimo-common.{KMER_LEN}mers.{NMUTS}mut.null")
config[section][f"fimo.{KMER_LEN}mer.{NMUTS}mut"] = JASPAR

crw.write(config, cfn)
#FA = os.path.join(OUTDIR, f"{CL}.{KMER_LEN}mers.{NMUTS}mut.nulls.fo.fa")

# run fimo
if os.path.exists(OUTDIR) is False:
    os.mkdir(OUTDIR)
    
os.chdir(OUTDIR)

cmd = [
    "fimo",
    JASPAR,
    FA, "&& mv" ,
    os.path.join(OUTDIR, "fimo_out"),
    os.path.join(OUTDIR, f"fimo-common.{KMER_LEN}mers.{NMUTS}mut.null")
    
]
print(" ".join(cmd))
os.system(" ".join(cmd))