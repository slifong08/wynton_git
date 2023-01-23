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

# In[ ]:


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


config_tag = sys.argv[1]
NBASES = sys.argv[2]

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

# nullomer GENCODE intersections
OVERLAP = config[f"DATAxGENCODE"]["OVERLAP"]  
NOOVERLAP_REF = config[f"DATAxGENCODE"]["NOOVERLAP_REF"]

# nullomer shuffles
EX_EXP_SHUF = config["SHUFFLE"]["ex_exp_star"]
NOEX_EXP_SHUF = config["SHUFFLE"]["noex_exp_star"]

# phylop files, wiggle track executables
PHYLOP_BW = config["PHYLOP"]["100WAY_BW"]


PHYLOP_BW = "/wynton/group/databases/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw"
PHYLOP_PATH = os.path.join(config["DATA"]["PATH"], "phylop")
SRC_BW = "/wynton/home/ahituv/fongsl/nullomers/src/bigWigSummary"

config["PHYLOP"]["100WAY_BW"] = PHYLOP_BW
config["PHYLOP"]["PATH"] = PHYLOP_PATH
config["PHYLOP"]["FLANK_BP"]=NBASES

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
            break

        for line in expanded:
            print("expanded", line)
            break
    
    return expanded, out


# ## parse bigwig for bed coordinates

# In[ ]:


def getBwBedLineSummary(bw, bed_line_list, nsteps, outfile, src):
    
    chr_, start_, end_, nid = bed_line_list
    
    ### insert a check for existing file, number of lines in file. 
    
    with open(outfile, "a") as writer:
        newline = "\t".join(bed_line_list)+"\t"
        writer.write(newline)
        writer.close()
                            
    cmd = f"{src} {bw} {chr_} {start_} {end_} {nsteps} >> {outfile}"
    
    os.system(cmd)
    


# ## expand + bw extract function

# In[ ]:


def expandAndBwExtract(bed, nbases, bw, outpath):
    
    # make the outfile
    name = (bed.split("/")[-1]).strip(".bed") + "-phylop.txt"
    outfile = os.path.join(outpath, name)

    # expand bed coordinates
    ex, outex = expand_bedcoor(bed, nbases)
    
    # get phylop values from bigwig
    for line in ex:
        bedlist = [line[0], line[1], line[2], line[3]]
        getBwBedLineSummary(bw, bedlist, nbases*2, outfile, SRC_BW)
        
    return outfile


def parallelBigWigVal(file_list, nbases, bw, outpath):
    
    """
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
    


# # Analyze 

# ## functions for analysis

# In[ ]:


def extract_phylop_values(bedfile):
    
    """
    return list of phylop average values from .bed file
    
    input
        bedfile (str) - path to bed file intersected w/ phylop.bw using bigWigAverageOverBed (above)
        
    method
        1. instantiate a list
        2. open the bedfile
        3. parse bedfile lines
        4. str split line to get relevant bedfile info, including phylop mean val at 4th index 
        5. append the value to a list of phylop values
    return
        phylop_vals (list) - list of mean phylop_values for all elements in the list. 
        
    """
    
    phylop_vals=[]  #1 collect phylop values from file
    
    #2 open the phylop.bed file
    with open(bedfile, "r") as phylofile:
        
        #3 parse lines
        for line in phylofile.readlines():
            
            #4 str split line to extract values
            chr_, posminusone, pos, nid, mean_phylop_val = (line.strip("\n")).split("\t")
            
            #5 append phylop_val to list of values
            phylop_vals.append(float(mean_phylop_val))
            
    return phylop_vals


# In[ ]:


def make_pdDataFrame(val_list, exp_val_list, sample_id):
    
    # make pd dataframe from val_list
    df=pd.DataFrame({"phylop":val_list})
    df["sample"] = sample_id
    
    # make pd dataframe from expected
    exp_df=pd.DataFrame({"phylop":exp_val_list})
    exp_df["sample"] = "shuffle"

    # concat the two dataframes together and reset index
    df=pd.concat([df, exp_df]).reset_index()
    
    return df


# In[ ]:


def results_stat_plot(df, sample_id):
    
    pp.fonts()
    
    # plot histplot
    fig,ax=plt.subplots(figsize=(6,6))
    sns.histplot(x="phylop", data=df, hue="sample", common_norm=False, stat="percent")
    ax.set(title=sample_id)
    
    # plot boxplot
    fig,ax=plt.subplots(figsize=(6,6))
    sns.boxplot(x="sample", y="phylop", data=df, notch=True, showfliers=False)
    ax.set(title=sample_id)

    # do MWU
    print(stats.mannwhitneyu(vals, exp_vals))

    # print descriptive stat
    print(df.groupby("sample")["phylop"].describe())




# # Main 

# ## exonic, non-exonic; expand; extract bw values

def main(argv):

    runlist =[
            OVERLAP,
            NOOVERLAP_REF
            ]

    for BED in runlist:

        print(BED)

        out = expandAndBwExtract(BED, NBASES, PHYLOP_BW, PHYLOP_PATH)

    # ## shuffled exonic, non-exonic; expand; extract bw values

    exp_lists = [
                glob.glob(EX_EXP_SHUF), 
                glob.glob(NOEX_EXP_SHUF)
                ]


    for exp_list in exp_lists:
        print(exp_list)
        exps_phylop = parallelBigWigVal(exp_list, NBASES, PHYLOP_BW, PHYLOP_PATH)
    

    # # write config

    crw.write_config(config, configname)

if __name__ == "__main__":
    main(sys.argv[1:])

    