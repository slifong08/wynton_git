#!/usr/bin/env python
# coding: utf-8

# 20221122
# sarahfong
# 
# 
# ### separates nullomers into coding and non-coding files
# - REQUIREMENT - GENCODE annotation for coding (user-specified)
# 
# - using GENCODE annotation bed 
# 
# - use pybedtools to intersect nullomers w/ GENCODE
# 
# ### separate non-coding into ref chromosomes and alternative haplotype chromosomes
# - drop redundant non-codingn coordinates (i.e. one locus, multiple nullomer options)
# 
# 
# ### preprocessing:
#     /wynton/home/ahituv/fongsl/nullomers/bin/pipeline/[0-3]*.py


import os
import pandas as pd
import pybedtools as pbt
import subprocess
import sys


# # config 

config_tag = sys.argv[1]


# append path
sys.path.append("/wynton/home/ahituv/fongsl/tools/py_/")

# import config reader
import config_readwrite as crw
import count_lines as cl

config_name = os.path.join(os.getcwd(), config_tag)

config, configfile_name = crw.read_config(config_name)


# In[5]:


# select config variables
ANNOT_GENCODE = config["GENCODE"]["ANNOT"]

GENCODE_BED = config["GENCODE"]["no_rmsk"]
NULLOMER_BED = config["DATA"]["no_rmsk"]


# # functions 

# In[7]:


def pbtIntersection(gencode_bed, test_bed, annot):
    
    out = test_bed.strip(".bed") + f".GENCODE-{annot}_overlap.bed"
    
    # make pybedtool objects
    gencode = pbt.BedTool(gencode_bed)
    test = pbt.BedTool(test_bed)
    
    # intersect -u to keep only unique loci from A
    
    test.intersect(gencode, u=True, output=out)
    

    print(cl.count_lines(out))
    print("\n\n", out)
    
    return out


# In[8]:


def pbtSubtraction(gencode_bed, test_bed, annot):
    
    out = test_bed.strip(".bed") + f".GENCODE-{annot}_no-overlap.bed"
    
    # make pybedtool objects
    gencode = pbt.BedTool(gencode_bed)
    test = pbt.BedTool(test_bed)
    
    # intersect -v to keep only the variants that DO NOT overlap (this is a subtraction. )
    # no -u for subtract. Need to assess unique values later
    
    test.intersect(gencode, v=True, output=out)
    
    print(cl.count_lines(out))
    print("\n\n", out)
    
    return out
    
    

def removeAltHap(bed_file):
    
    """
    return .bed file with only reference chromosome nullomers
        ***specifically for the non-coding file, coding file (from GENCODE) is ref chromosome only. 
    
    input
        bed_file (str) - path to bed file
        
    method
        1. if ref_bed output file does not exist, make one. 
        
        2. open non-coding file
        3. make list of reference chromosome names
        4. subset non-coding dataframe for reference chromosomes. 
        5. write REF_BED

    return 
        REF_BED (str) - bed file w/ only reference chromosomes.

        
    """
    

    REF_BED = bed_file.strip(".bed") + "-ref.only.bed"
    
    if os.path.exists(REF_BED) is False:
        # open the non-coding file (gencode overlap subtracted mutations)

        df = pd.read_csv(bed_file, sep='\t', header=None)
        print()

        # make a list of the reference chromosomes
        from chr_functions import make_chr_list

        ref_chr_list = make_chr_list()
        ref_chr_list.extend(["chrX", "chrY"])


        # ## subset non-coding reference chromosomes

        # subset reference chromosome dataframe
        ref = df.loc[df[0].isin(ref_chr_list)]
        
        # drop dups
        ref = ref[[0,1,2,3]].drop_duplicates()
        print("before", df.shape[0], "after", ref.shape[0], "diff", df.shape[0]-ref.shape[0])

        # how many mutations are on alternative haplotypes? 

        # How many nullomers are at alternative haplotype loci?
        # 
        #     1866/11317 mutations are in alternative haplotypes (17%)(all GENCODE)
        #     1868/24632 mutations are in alternative haplotypes (7.6%) (exons only)
        # 
        #     ... the rest of the mutations are NON-CODING and in REFERENCE chromosome coordinates. 


        # save reference chromosome non-coding mutation dataframe
        ref.to_csv(REF_BED, sep='\t', header=False, index=False)
    
    return REF_BED


    
def main(argv):


    # # intersect mutations w/ GENCODE 

    # ### write GENCODE_OVERLAP


    # intersection
    GENCODE_OVERLAP = pbtIntersection(GENCODE_BED, NULLOMER_BED, ANNOT_GENCODE) # remove duplicates


    # ### write GENCODE_NOOVERLAP


    # subtraction
    GENCODE_NOOVERLAP = pbtSubtraction(GENCODE_BED, NULLOMER_BED, ANNOT_GENCODE)


    # remove alternative haplotypes
    GENCODE_NOOVERLAP_REF = removeAltHap(GENCODE_NOOVERLAP)

    crw.check_section(config, "DATAxGENCODE")
    config[f"DATAxGENCODE"]["OVERLAP"] = GENCODE_OVERLAP   # write
    config[f"DATAxGENCODE"]["NOOVERLAP"] = GENCODE_NOOVERLAP   # write
    config[f"DATAxGENCODE"]["NOOVERLAP_REF"] = GENCODE_NOOVERLAP_REF   # write

    # write the config file
    crw.write_config(config, configfile_name)

    print(configfile_name)
    
if __name__ == "__main__":
    main(sys.argv[1:])





