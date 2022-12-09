#!/usr/bin/env python
# coding: utf-8

# 20221122
# sarahfong
# 
# ### Format mutation file
# 
# add 
# - nullomer_iD
# 
# 
# write:
# - vcf-like file (#chrom, pos, id, ref, alt
# - bed-like file (#chr, pos-1, pos, ref, alt)

# In[1]:


import os
import pandas as pd
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


# In[3]:


# select config variables

MUTATION = config["DATA"]["MUTATIONS"]


# # functions 

# In[4]:


def makeNulloID(df):
    
    """
    add unique nullomer id tag to each uniq nullomer locus
    
    input
        df (pd.DataFrame object) of nullomer mutation list
    
    output 
        nullo_df (od.DataFrame object) of uniq nullomer loci w/ 
    """
    
    nullo_df = df[[0,1]].drop_duplicates().sort_values(by=[0,1]) # get unique loci [CHR, POS], sort 

    print(df.shape[0], nullo_df.shape[0])
    
    nullo_df["NID"] = "N_" +nullo_df.index.map(str) # add NULLOMER ID
    
    return nullo_df


def main(argv):
    
    # files to write
    VCF = MUTATION.strip(".txt") + ".vcf"
    BED = MUTATION.strip(".txt") + ".bed"
    BED_UNQ = MUTATION.strip(".txt") + "-uniq.bed"
    
    # # mutations to vcf, bed 

    df = pd.read_csv(MUTATION, sep='\t', header=None)

    # how big is the dataframe? 
    print(df.shape)

    nullo_id_df = makeNulloID(df)  # add nullomer ids based on uniq chrom/pos
    df = pd.merge(df, nullo_id_df, how = "left")

    df[5] = df[1]-1 # POS-1, or START position
    df[6] = "chr"+df[0]
    df.head()


    # ### make the VCF

    vcf = df[[0,1,"NID",2,3]].sort_values(by=[0,1])

    # save vcf-like file
    vcf.to_csv(VCF, sep='\t', index=False, header=False)
    vcf.head()


    # ### make bed


    bed = df[[6,5,1,"NID",2,3]].sort_values(by=[6,1])

    # save bed-file
    bed.to_csv(BED, sep='\t', index=False, header=False)
    print(bed.shape)
    bed.head()


    # ## uniq loci

    bed_uniq = df[[6,5,1,"NID"]].sort_values(by=[6,1]).drop_duplicates()

    # save bed-file
    bed_uniq.to_csv(BED_UNQ, sep='\t', index=False, header=False)
    print(bed_uniq.shape)
    bed_uniq.head()

    # write to config
    config["DATA"]["VCF"] = VCF  # write
    config["DATA"]["BED"] = BED  # write
    config["DATA"]["BED_UNQ"] = BED_UNQ  # write


    # write the config file
    crw.write_config(config, configfile_name)

    print(configfile_name)
    
if __name__ == "__main__":
    main(sys.argv[1:])



