#!/usr/bin/env python
"""
make non-coding exlucsion background from blacklist+gencode+
"""
import glob
import os
import subprocess
import sys


# # config 
config_tag = sys.argv[1]  # str
ITERS = sys.argv[2]  #str


# append path
sys.path.append("/wynton/home/ahituv/fongsl/tools/py_/")

# # config 

# import config reader
import config_readwrite as crw

config_name = os.path.join(os.getcwd(), config_tag)

config, configfile_name = crw.read_config(config_name)

# select config variables

BUILD = config["DATA"]["BUILD"]
DATA_PATH = config["DATA"]["PATH"]

ANNOT = config["GENCODE"]["ANNOT"]
GENCODE = config["GENCODE"]["BED"]
FLAT_GENCODE = config["GENCODE"]["MERGED"]

RMSK_BED = config["RMSK"]["BED"]

# files to shuffle. 
GENCODE_OVERLAP = config["DATAxGENCODE"]["OVERLAP"]
GENCODE_NOOVERLAP = config["DATAxGENCODE"]["NOOVERLAP_REF"]


PY = "/wynton/home/ahituv/fongsl/tools/genome/shuf_wynton-w_config.py"
QSUB = "/wynton/home/ahituv/fongsl/tools/qsub/shuf_nullomer-w_config.sh"

# # functions 


def makeBkgdExclusion(data_path, build, gencode, rmsk):
    
    """
    make .bed combining blacklist, GENCODE coding regions, repeatmasker
    
    input
    
    method
       
    return
        out (str) - path to background exclusion .bed file
        blacklist (str) - blacklist file. 
        
    """
    
    BLACK_LIST = f"/wynton/home/ahituv/fongsl/dna/black_list/{build}-blacklist.v2.bed"
    OUT = os.path.join(data_path, f"bkgdExclusions.blacklist-{build}.rmsk.gencode.bed")
    
    if os.path.exists(OUT) is False:
        print(" ".join(["cat", BLACK_LIST, rmsk, gencode, "| cut -f 1,2,3 >",  OUT]))
        subprocess.run(["cat", BLACK_LIST, rmsk, gencode, "| cut -f 1,2,3 >",  OUT], stdout=subprocess.PIPE)
    
    return OUT, BLACK_LIST

def writeShuffleConfig(bed, iters, build, outpath, incl, section, annotation):
    
    """
    return config name, file to write to 
    
    inputs
        bed (str) - input bed
        iters (int) - n iters
        build (str) - genome build
        outpath (str) - path to deposit shuffles
        incl (str) - path to file to include/exclude in shuffles. 
        
    method
        1. get path
        2. get config name, file name
        3. touch config file
    
    return 
        configname (str) - "config-<<type>>.ini"
        configfile (str) - FULL PATH + "config-<<type>>.ini"
        
    """
    #1
    path = os.getcwd()
    
    #2
    configname = "config-shuf-" + annotation +".ini"
    configfile = os.path.join(path, configname)
    
    #3
    cmd = "touch "+ configfile
    if os.path.exists(configfile) is False:
        os.system(cmd)
        
    config, configfile_name = crw.read_config(configfile)
    
    crw.check_section(config,section)
    
    config[section]["BED"] = bed
    config[section]["ITERS"] = iters
    config[section]["BUILD"] = build
    config[section]["OUTDIR"] = outpath
    config[section]["INCL"] = incl

    # write the config file
    crw.write_config(config, configfile_name)
    
    return configfile_name

def shuffle(config, section, cmd):
    
        
    """
    launch qsub command to run shuffles. 
    
    input
        mutation_bed (str) - bed file to shuffle
        background (str) - path to file to include or exclude
        iters (str) - number of shuffles
        shuf_path (str) - path to shuffles 
    
    method
        run qsub as subprocess. 
       
    return
        out (str) - path to background exclusion .bed file
        blacklist (str) - blacklist file. 
        
    """
    #### STOPPED HERE. NEED TO SPLIT FILE NAMES AND CREATE WYNTON QSUB COMMAND TO SHUFFLE. 
    if ".py" in cmd:
        print(" ".join(["python3", cmd, config, section]))
        subprocess.run(["python3", cmd, config, section], stdout=subprocess.PIPE)
    else:
        print(" ".join(["qsub", cmd, config, section]))
        subprocess.run(["qsub", cmd, config, section], stdout=subprocess.PIPE)

def main(argv):
    
    # get RMSK.bed, write noncoding background exclusion file. 
    EXCL_BKGD, BLACK_LIST = makeBkgdExclusion(DATA_PATH, BUILD, GENCODE, RMSK_BED)
    
    # sections, files to write
    sections = ["BKGD", "SHUFFLE"]

    for s in sections:
        crw.check_section(config, s)
                       
    # write to config
    config["BKGD"]["EXCL"] = EXCL_BKGD  # write
    config["BKGD"]["BLACKLIST"] = BLACK_LIST  # write


    # shuffle the datasets
    SHUFFLE_PATH = os.path.join(DATA_PATH + "shuffle")

    # write the shuffle config

    shuf_config = writeShuffleConfig(GENCODE_OVERLAP, ITERS, BUILD, SHUFFLE_PATH, FLAT_GENCODE, "OVERLAP", ANNOT)
    shuf_config = writeShuffleConfig(GENCODE_NOOVERLAP, ITERS, BUILD, SHUFFLE_PATH, EXCL_BKGD, "NO-OVERLAP", ANNOT)
    
    # write to the config
    
    config["SHUFFLE"]["PATH"] = SHUFFLE_PATH
    config["SHUFFLE"]["ITERS"] = str(ITERS)
    config["SHUFFLE"]["QSUB"] = QSUB
    config["SHUFFLE"]["PY"] = PY
    config["SHUFFLE"]["EX_EXP_STAR"] = os.path.join(SHUFFLE_PATH, "shuf-" + (GENCODE_OVERLAP.split("/")[-1]).strip(".bed") +"-*.bed") 
    config["SHUFFLE"]["NOEX_EXP_STAR"] = os.path.join(SHUFFLE_PATH, "shuf-" + (GENCODE_NOOVERLAP.split("/")[-1]).strip(".bed") +"-*.bed")     
    
    # run shuffles
    shuffle(shuf_config, "OVERLAP", PY)
    shuffle(shuf_config, "NO-OVERLAP", PY)
    
    # write the config file
    crw.write_config(config, configfile_name)

    print(configfile_name)
    
if __name__ == "__main__":
    main(sys.argv[1:])