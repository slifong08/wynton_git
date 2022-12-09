#!/usr/bin/env python
# coding: utf-8

# 20221123
# 
# sarahfong
# 
# ### intersect nullomers with phastCons
# 
# phastCons100way hg38 .bed file was downloaded from the UCSC genome table browser
# 
# 
# 1. count the number of overlapping nullomers. 
#     Remove duplicate nullomers (i.e. where one position can have two nullomers)
# 2. calculate expectation and intersect w/ phastCons
# 
# ### preprocessing
#     /wynton/home/ahituv/fongsl/nullomers/bin/GENCODE/0_download_format_gencode.ipynb
# 
#     /wynton/home/ahituv/fongsl/nullomers/bin/GENCODE/0_format_mutation_file.ipynb
# 
#     /wynton/home/ahituv/fongsl/nullomers/bin/GENCODE/1_separate_coding_non-coding.ipynb

# In[1]:


import glob
from joblib import Parallel, delayed
import os
import pybedtools as pbt
import subprocess
import sys

import numpy as np
import seaborn as sns
import datetime


config_tag = sys.argv[1]


# append path
sys.path.append("/wynton/home/ahituv/fongsl/tools/py_/")

# import config reader
import config_readwrite as crw
import count_lines as cl

config_name = os.path.join(os.getcwd(), config_tag)

config, configname = crw.read_config(config_name)


# In[3]:


# select config variables
ANNOT = config["GENCODE"]["ANNOT"]
GENCODE = config["GENCODE"][f"{ANNOT}_BED"]
FLAT_GENCODE = config["GENCODE"]["MERGED"]

# nullomer intersections
OVERLAP = config[f"DATAxGENCODE"]["OVERLAP"]  
NOOVERLAP_REF = config[f"DATAxGENCODE"]["NOOVERLAP_REF"]

# shuffles
DATA_PATH = config["DATA"]["PATH"]
SHUF_PATH = config["SHUFFLE"]["PATH"]
RESULTS_PATH = config["RESULTS"]["PATH"]

#PHASTCONS = config["PHASTCONS"]["100WAY"] # write
#EX_EXP = config["PHASTCONS"]["EXON_EXP"]  # write
#NOEX_EXP = config["PHASTCONS"]["NOEXON_EXP"]  # write


PHASTCONS = "/scratch/fongsl/nullomers/data/phastCons/phastCons100way_hg38.bed"


def phastcons_intersection(phastcons_bed, test_bed):
    """
    intersect test bed w/ phastcons elements
    
    input
        phastcons_bed (str) - path to phastcons.bed file
        test_bed (str) - path to test.bed file
        
    method
        1. turn bed files into pybedtool objects
        2. intersect bed and phastcons files
        3. count number of overlaps w phastcons
        
    return
        test_int (pybedtool object) - intersected pybedtools object. 
        count (int) - count of lines in intersection
        
    """
    
    #1
    phast = pbt.BedTool(phastcons_bed)
    test = pbt.BedTool(test_bed)
    
    #2
    test_int = test.intersect(phast)
    
    count = sum(1 for line in test_int) 
    print(count)
    
    #3
    return count   # return pbt.object, count



def parallel_intersections(shuffle_list, phastcons):

    #num_cores = multiprocessing.cpu_count()
    num_cores = 16
    print("number of cores", num_cores)

    # run parallel jobs

    exp = Parallel(n_jobs=num_cores, verbose=100, prefer="threads")(delayed(phastcons_intersection)(phastcons, shuf_iter) for shuf_iter in shuffle_list)
    
    return exp


# In[10]:


def exp(bed, incl, annot, phastcons, shuf_path):
    
    
    if incl is None:
        
        shufs = glob.glob(os.path.join(shuf_path, f"shuf*-{annot}_no-overlap*.bed"))
    else:
        shufs = glob.glob(os.path.join(shuf_path,f"shuf*-{annot}_overlap-*.bed"))
    print( "n shuffles to intersect", len(shufs))
        
    # parallel process here.    
    exp = parallel_intersections(shufs, phastcons)
   
    return exp


# In[11]:


def write_expectation(outfile, exp_list):

    with open(outfile, "w") as results:
        
        for i in exp_list:
            line = f"{i}\n"
            results.write(line)
        
        results.close()


        
def calculateEmpiricalP(obs, exp_sum_list):
    """
    return two lists
        (1) info - vector w/  
                n_obs, 
                median_exp, 
                std, 
                fold-change  # calculated from the median of expected shuffle 
                p_val
                
        (2) fold_changes- vector expected fold changes (to calculate confidence interval)
        
    input
        observed overlap count (int)
        list of expected overlap counts (list of ints)
    
    method
        1. get median of expected overlap counts
        2. get standard deviation of expected overlap counts
        3. center expected overlap counts at median
        4. Sum the number of centered expected counts greater than observed centered count
            This is two tailed because it evaluates both sides of the distribution (w/ abs value). 
        5. calculate fold change as observed/ median expected w/ pseudo count
        6. calculate fold change of each "obs"/ expected w/ pseudo count
        7. calculate the p-value as count of equal or more extreme values than observed value
        8. return list of empirical info + fold changes
        
        
    
    """
    #1
    mu = np.median(exp_sum_list)  # median of exp.dist
    
    #2
    sigma = np.std(exp_sum_list)  # std
    
    #3
    dist_from_mu = [exp - mu for exp in exp_sum_list] # center the distribution 
    
    #4
    p_sum = sum(1 for exp_dist in dist_from_mu if abs(exp_dist) >= abs(obs - mu)) # count values >= centered obs

    #5
    fold_change = (obs + 1.0) / (mu + 1.0) # fold change obs from median expected w pseudo count
    
    #6
    fold_changes = list((obs + 1.0) / (m + 1.0) for m in exp_sum_list) # fold change obs from /each exp w pseudo count
    
    #7
    p_val = (p_sum + 1.0) / (len(exp_sum_list) + 1.0)  # probability of observing obs-like value equal or more extreme in expected distribution
    
    #8
    info = [
            obs, 
            mu, 
            sigma, 
            fold_change, 
            p_val, 
            str(datetime.datetime.now())
            ]
    
    return info, fold_changes

def writeStats(outfile, stat_list_tup):
    """
    append stat to stat file. 
    
    input 
        outfile (str) - full path to outfile.txt
        stat_list (tuple) - (description of stat, stat_value) 
        
    method
        1. if outfile (or directory) does not exist, make outfile
        2. write tuple line into outfile. 
        3. close file
        
    """
    #1
    if os.path.exists(outfile) is False:
        
        if os.path.dirname(outfile) is False:
            os.system(f"mkdir {os.path.dirname(outfile)}")
        
        os.system(f"touch {outfile}")
    #2    
    with open(outfile, "a") as results:
        
        for descriptor, stat in stat_list_tup:
            line = f"{descriptor}\t{stat}\n"
            results.write(line)
        #3
        results.close()
    

def main(argv):

    # variables for writing results
    RE = os.path.join(RESULTS_PATH, "phastcons")
    RE_STATS = os.path.join(RE, "phastCons.stats.txt")
    EX_EXP = os.path.join(DATA_PATH, f"{ANNOT}_exp-counts.txt")
    NOEX_EXP = os.path.join(DATA_PATH, f"no-{ANNOT}_exp-counts.txt")

    # check config for PhastCons
    crw.check_section(config, "PHASTCONS")
    
    config["PHASTCONS"]["100WAY"] = PHASTCONS
    config["PHASTCONS"][f"{ANNOT}_EXP"] = EX_EXP # write
    config["PHASTCONS"][f"NO-{ANNOT}_EXP"] = NOEX_EXP # write
    config["RESULTS"]["PHASTCONS"] = RE
    config["RESULTS"]["PHASTCONS_STATS"] = RE_STATS
    results_list = []

    # ## Exonic x phastcons


    # coding x phastcons intersection    
    obs_coding = phastcons_intersection(PHASTCONS, OVERLAP)
    n_total_coding=(sum(1 for line in open(OVERLAP, "r")))
    
    results_list = []
    results_list.extend([(f"n {ANNOT} total", n_total_coding), 
                        (f"n {ANNOT} w/ phastcons overlap", obs_coding), 
                        (f"% {ANNOT} w/phastcons overlap", round(obs_coding/n_total_coding, 2))
                        ])
    
    writeStats(RE_STATS, results_list)

    # How many EXON nullomer loci overlap phastcons element? 
    # 
    #     2150/4593 exon loci are conserved in phastCons (46.8%) (ANY GENCODE EXON)
    #     
    #     OLD - 2891/19807 coding loci are conserved in phastCons (15%) (ANY GENCODE)

    # # empirical expectation

    # ## exonic
    # - shuffle in coding background (GENCODE)
    # - FIRST, need to flatten GENCODE coordinates. 


    # # expectation 

    # ## exonic

    # ### shuffle exonic nullomers into GENCODE background

    exp_Overlap = exp(OVERLAP, FLAT_GENCODE, ANNOT, PHASTCONS, SHUF_PATH)
    write_expectation(EX_EXP, exp_Overlap)


    # In[46]:

    results = calculateEmpiricalP(obs_coding, exp_Overlap)

    results_header = [
                    "obs", 
                    "mu", 
                    "sigma", 
                    "fold_change", 
                    "p_val", 
                    "run",
                    ]
    
    # write fold change for coding
    writeStats(RE_STATS, zip(results_header, results))

    # ## non-coding x phastcons
    

    # non-coding x phastcons intersection    
    obs_noncoding = phastcons_intersection(PHASTCONS, NOOVERLAP_REF)

    n_total_noncoding = sum(1 for line in open(NOOVERLAP_REF, "r"))


    # How many NON-EXON nullomer loci overlap phastcons element? 
    #     
    #     1178/24546 non-coding loci are conserved (4.8%) (ANY GENCODE EXON)
    #     
    #     OLD - 439/9423 non-coding loci are conserved (5%) (ANY GENCODE)

    

    results_list = []
    results_list.extend([(f"n non-{ANNOT} total", n_total_noncoding), 
                        (f"n non-{ANNOT} w/ phastcons overlap", obs_noncoding), 
                        (f"% non-{ANNOT} w/phastcons overlap", round(obs_noncoding/n_total_noncoding, 2))
                        ])
    print(results_list)
    writeStats(RE_STATS, results_list)


    exp_noOverlap = exp(NOOVERLAP_REF, None, ANNOT, PHASTCONS, SHUF_PATH)
    write_expectation(NOEX_EXP, exp_noOverlap)


    results = calculateEmpiricalP(obs_noncoding, exp_noOverlap)
    
    writeStats(RE_STATS, zip(results_header, results))
    
    
if __name__ == "__main__":
    main(sys.argv[1:])