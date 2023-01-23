#!/usr/bin/env python
# coding: utf-8

# 20221129
# 
# sarahfong
# 
# ### intersect nullomers, empirical shuffle with phylop 100way bigWig
# 
# split by exonic/non-exonic
# 
# 
# use bigWigSummary executable from UCSC to get phylop 
# 
# 
# compare nullomers v. empirical background

import glob
from joblib import Parallel, delayed
import os
import pybedtools as pbt
import subprocess
import sys


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import statsmodels as sm

# append path
sys.path.append("/wynton/home/ahituv/fongsl/tools/py_/")

# import config reader
import config_readwrite as crw
import zippery

config_tag = sys.argv[1]
NBASES = sys.argv[2]
BED = sys.argv[3]

print(BED, type(BED))


# append path
sys.path.append("/wynton/home/ahituv/fongsl/tools/py_/")

# import config reader
import config_readwrite as crw
import count_lines as cl
import plot_params as pp

config_name = os.path.join(os.getcwd(), config_tag)

config, configname = crw.read_config(config_name)


# select config variables
ANNOT = config["GENCODE"]["ANNOT"]  
DATA_PATH = config["DATA"]["PATH"]


# phylop files, wiggle track executables
PHYLOP_BW = config["PHYLOP"]["100WAY_BW"]

PHYLOP_BW = "/wynton/group/databases/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw"
PHYLOP_PATH = os.path.join(DATA_PATH, "phylop")
SRC_BW = "/wynton/home/ahituv/fongsl/nullomers/src/bigWigSummary"

config["PHYLOP"]["100WAY_BW"] = PHYLOP_BW
config["PHYLOP"]["PATH"] = PHYLOP_PATH
config["PHYLOP"]["FLANK_BP"]= str(NBASES)

crw.check_section(config, "SRC")
config["SRC"]["bigwigsummary"] = SRC_BW


# # functions


def expand_bedcoor(bedfile, nbases):
    """
    write expanded bedcoordinates file. 
    
    inputs
        bedfile (str) - full path to bed file to expand
        nbases (int) - number of bases to expand by IN EACH DIRECTION
        
    method
        
        1. build outfile str
        2. make pybedtool object (if outfile does not exist, else skip)
        3. use bedtools slop to expand bed coordinates, save output
        4. unit test expansion
        
    returns
        expanded (bedtool object) - expanded bed coordinates
        outbed (str) - full path to expanded coordinate bedfile
    """
    
    #1
    out = bedfile.strip(".bed") + f"-EXPANDED-{nbases}.bed"
    
    if os.path.exists(out) is False:
        #2
        bed = pbt.BedTool(bedfile)

        #3
        expanded = bed.slop(genome="hg38", b=nbases, output=out)
    else:
        print("expanded", nbases, " in either direction, already")
        expanded = pbt.BedTool(out)
    
    #4
    with open(bedfile, "r") as file:
        for line in file:
            print("original", line)
            break  # only print the first original line in the file

        for line in expanded:
            print("expanded", line) 
            break  # only print he first expanded line in the file
    
    return expanded, out


# ## parse bigwig for bed coordinates



def getBwBedLineSummary(bw, bed_line_list, nsteps, outfile, src):
    
    """
    write bed coordinates + phylop values for a bed-coordinate region. 
    
    input
        bw (str) - path to bigwig file to parse
        bed_line_list (list) - list with CHR, START, STOP, and 4th column annotations. 
            Should correspond to 1 .bed line.
        nsteps (int) - steps to summarize the phylop sequence in. 
            e.g. If 1000 bases in sequences and nsteps == 1000, will return each bp phylop score. 
            e.g. if 1000 bases in sequences and nsteps == 10, will average (?) 100bp phylop scores into 1 output. 
        outfile (str) - full path to write file. 
        src (str) - path to ./bigWigSummary executable. 
        
    method
        1. unpack bed_line_list into CHR, START, STOP, and sample ID variables
        2. run bigWigSummary command, get phyloP scores for the bed coordinates (output). 
        3. make new line (str) w/ bedcoordinates, id, phylop results
        4. write bed coordinates and phylop to the outfile
        
        
    """
    
    #1
    chr_, start_, end_, nid = bed_line_list
    
    # TO DO: insert a check for existing file, number of lines in file. Concerned about parallelized runs and writing to wrong coordinates. 
    
    #2              
    cmd = f"{src} {bw} {chr_} {str(start_)} {str(end_)} {str(nsteps)}"
 
    results = subprocess.getoutput(cmd)
    
    #3 make a str (tab-delimited) of the bed coordinates and phylop results
    newline = "\t".join(bed_line_list) +"\t"+ results.replace(" ", "\t")+"\n"
    
    #4 
    with open(outfile, "a") as writer:
        writer.write(newline)
        writer.close()


# ## expand + bw extract function

def getSummaryStats(phylop_file):
    
    """
    return summary stats for each position in a phylop window
    
    input
        phylop_file (str) - path to phylop bw extracted file from above run
        
    method
        1. make a pandas dataframe out of the dictionary
            
            col 0-3 = CHR, START, STOP, Nullomer ID
            col 4-1003 = each bp position phylop score flanking nullomer
            
            col 503 = nullomer. 

        2. filter
            2.1 remove any lines w/ no data 
            2.2 replace "n/a" str value
        
        3. turn dataframe into long form data, where 1 column is the position, and another column is the phylop value
            index on CHR, START, STOP, Nullomer ID
        
        4. compute summary stats per position across nullomer sequence for phylop scores.
            4.1 also correct for
        
        5. write summary stats to outfile
    
    return 
    
        summary stats outfile (str) - write summary stats as file, 
            ## so that stats can be done on empirical summary statistics. 
            
        
    """
    out = phylop_file.strip(".txt") + "-SUMMARY_STATS.txt"

    #1 make a pandas dataframe
    df = pd.read_csv(phylop_file, sep = '\t', header = None, low_memory = False)
    

    #2
    df = df.loc[df[4] != "no"] ## there are some cases where there is no data. 
    df = df.replace("n/a", None)
    print(df.shape)
    print(df.head())

    
    #3
    melted = pd.melt(df, id_vars= df.columns[:4], var_name="pos", value_name="phylop")
    print(melted.head())
    
    melted["phylop"] = melted["phylop"].astype(float)  # change data type
    
    #4
    summary_stats = melted.groupby("pos")["phylop"].describe().reset_index()
    
    summary_stats["pos"] = summary_stats["pos"] - 4 # correct bp position 
    
    #5 
    summary_stats.to_csv(out, sep ='\t', index=False)
    
    zippery.rezip_file(phylop_file)
    
    return out


def expandAndBwExtract(bed, nbases, bw, outpath):
    """
    combined function that:
        0. makes an outfile. 
            0.1 if scores have not been extracted already
            0.2 if scores have been and the file is zipped
        1. takes locus and expands the bed coordinates of the locus. 
        2. extracts the phylop base scores for the expanded locus, writing the scores as a vector
        3. write summary stats of file. (Summary stats at each bp)
            

    return 
        outfile (str) - name to output values. 
    """
    #0 make the outfile
    name = (bed.split("/")[-1]).strip(".bed") + "-phylop.txt"
    outfile = os.path.join(outpath, name)
    zipped_file = outfile + ".gz"

    #0.1
    if os.path.exists(outfile) is False and os.path.exists(zipped_file) is False:

        #1 expand bed coordinates
        ex, outex = expand_bedcoor(bed, nbases)
    
        #2 get phylop values from bigwig and write to the outfile. 
        for line in ex:
            bedlist = [line[0], line[1], line[2], line[3]]
            getBwBedLineSummary(bw, bedlist, nbases*2, outfile, SRC_BW)
    #0.2
    elif os.path.exists(zipped_file) is True:
        outfile=zippery.unzip_file(zipped_file)

    # 3
    summary_outfile = getSummaryStats(outfile)
    
    return summary_outfile


def parallelBigWigVal(file_list, nbases, bw, outpath):
    
    """
    OLD. FOR SHUFFLES ONLY.
    run parallel jobs to parse bigwigs for bed values

    input
        file_list (list) - list of strs w/ full paths to bed files to run in parallel
        nbases (int) - number of bases around nullomer mutation to get bw values from
        bw (str) - full path to bigwig file to parse
        outpath (str) - path to write outfile

        
    method
        run parallel jobs to extract plylop values from bw
        
    return
        exps_phylop (list) - list of outfiles
    
    notes
        ncores is fixed at 16
        
    """
    ncores=16
    
    exps_phylop = Parallel(
                        n_jobs=ncores, verbose=100, prefer="threads")\
                        (delayed(expandAndBwExtract)\
                        (bed, nbases, bw, outpath) for bed in file_list)
    return exps_phylop

# # Main 

# ## exonic, non-exonic; expand; extract bw values

def main(argv):

    
    # combined function that expands mutations and extracts phylop scores
    out = expandAndBwExtract(BED, NBASES, PHYLOP_BW, PHYLOP_PATH)
    

    # # write config

    crw.write_config(config, configname)

if __name__ == "__main__":
    main(sys.argv[1:])

        