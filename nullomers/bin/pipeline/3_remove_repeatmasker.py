#!/usr/bin/env python
# coding: utf-8

# 20221206
# sarahfong
# 
# REMOVE REPEATMASKER FROM:
#  Nullomer file
#  GENCODE file
# 
# ADD REPEATMASKER TO:
#   Non-coding blacklist background file. 
# write:
# - gencode bedfile w/ subtracted rmsk
# - bed file w/ subtracted rmsk
# - genomic background file w/ added rmsk

# In[1]:


import os
import pandas as pd
import pybedtools as pbt
import subprocess
import sys


# # config 
config_tag = sys.argv[1]

# append path
sys.path.append("/wynton/home/ahituv/fongsl/tools/py_/")

# # config 

# import config reader
import config_readwrite as crw
import zippery

config_name = os.path.join(os.getcwd(), config_tag)

config, configfile_name = crw.read_config(config_name)


# select config variables

DATA_PATH = config["DATA"]["PATH"]
MUTATION = config["DATA"]["BED_UNQ"]
GENCODE = config["GENCODE"]["MERGED"]
BUILD = config["DATA"]["BUILD"]


print("BUILD", BUILD)
# # functions 

# In[4]:


def format_repeatmasker(data_path, build):
    
    """
    make .bed from current repeatmasker annotations
    
    input
       data_path (str) - full path to where to write repeatmasker.bed file
       build (str) - genome build
    
    method
        run zmore on wynton rmsk file and save to a bed file. 
            6,7,8 = chr, start, stop,
            10 = strand
            11,12,13 = TE species, subfamily, family

    return
        rmsk_bed (str) - bed file written from wynton_rmsk
        wynton_rmsk (str) - original wynton rmsk file. 
        
    """
    
    WYNTON_RMSK = f"/wynton/group/databases/goldenPath/{build}/database/rmsk.txt.gz"
    LOCAL_RMSK = os.path.join(data_path, f"rmsk.txt.gz") 
    RMSK_BED = os.path.join(data_path, f"rmsk-{build}.bed")
    
    if os.path.exists(RMSK_BED) is False:
        mv_cmd = f"cp {WYNTON_RMSK} {data_path}"
        subprocess.call(mv_cmd, shell=True)
        
        unzipped = zippery.unzip_file(LOCAL_RMSK)
        cmd = f"cut -f 6,7,8,10,11,12,13 {unzipped} > {RMSK_BED}"

        subprocess.call(cmd, shell=True)
        
        zippery.rezip_file(unzipped)

    return RMSK_BED, WYNTON_RMSK


def subtractRmsk(bedfile, rmsk_bed):
    """
    return bedfile wout repeatmasker overlaps
    
    input
        bedfile (str) - path to bed file to subtract rmsk from
        rmsk_bed (str) - path to rmsk.bed file
        
    method
        1. make subtraction file. 
        2. if subatraction file does not exist 
            2.1 use pybedtools to subtract rmsk.bed
    
    return 
        subtraction (str) - path to subtraction file.  
    """
    
    #1 make the outfile
    SUBTRACTION = bedfile.strip(".bed") + "-woRMSK.bed"
    
    #2 if the file does not exist
    if os.path.exists(SUBTRACTION) is False:
        
        #2.1
        bed = pbt.BedTool(bedfile)
        rmsk = pbt.BedTool(rmsk_bed) 
        
        bed.subtract(rmsk, output=SUBTRACTION)
    
    return SUBTRACTION
        

def main(argv):
    
    # get RMSK.bed 
    RMSK_BED, WYNTON_RMSK = format_repeatmasker(DATA_PATH, BUILD)

    # subtract RMSK.bed from bedfiles:
    MUTATION_NORMSK = subtractRmsk(MUTATION, RMSK_BED)
    
    GENCODE_NORMSK = subtractRmsk(GENCODE, RMSK_BED)
    
    # files to write
    
    # write to config
    if config.has_section("RMSK") is False:
        config.add_section("RMSK") 
        
    config["RMSK"]["WYNTON"] = WYNTON_RMSK  # write
    config["RMSK"]["BED"] = RMSK_BED  # write
    config["DATA"]["NO_RMSK"] = MUTATION_NORMSK  # write
    config["GENCODE"]["NO_RMSK"] = GENCODE_NORMSK  # write


    # write the config file
    crw.write_config(config, configfile_name)

    print(configfile_name)
    
if __name__ == "__main__":
    main(sys.argv[1:])



