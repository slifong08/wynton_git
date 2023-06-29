#!/usr/bin/env python
# coding: utf-8

# In[2]:


from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import reverse_complement, Seq
from Bio.SeqUtils import gc_fraction
import config_readwrite as crw

import glob
import gzip
import numpy as np
import os, sys
import pandas as pd
import subprocess as sp
import time


# # find coordinates of kmers

# In[3]:


CL, MER, NMUTS, BUILD, FO_ONLY = "common", "15", "2", "hg38", True
config, cfn = crw.read(os.path.join(os.getcwd(), "config.ini"))


# In[4]:


# read

# concatenated nullomers
CONCAT= config[f"{CL}.{MER}mer.{NMUTS}mut"]["concat-13"]

# randomly sampled nullomers for single inserts
NULLS = config[f"{CL}.{MER}mer.{NMUTS}mut"]['null_seqid_greaterthan1bp_sample']

# bed file of common MPRA coordinates
BED = config[CL][f"bed_{BUILD}"]

# extended 4096bp fasta file for common MPRA coordinates
FA_EXT_BED = config[CL][f"fa_extended_{BUILD}"]

# extended 4096 bed file for common MPRA coordinates
EXTENDED = config[CL][f"bed_extended_{BUILD}"]

# fast file for common MPRA coordinates
FA = config[CL][f"fasta_{BUILD}"]

# scaffolds for inserting single nullomers. 
# Chosen because active (top 97.5% Agarwal) | inactive (IQR Agarwal) across 3 cell lines.
# Randomly sample 150 elements from these two categories. 
# Predict scaffold activity with Sei
# Pick scaffolds w/ predicted activity in Sei. (Are these sequences more active than not predicted?)
JOINT_LIBRARY=config["MPRA_AGARWAL"]["joint_library"]
SCAFFOLDS = config["MPRA_AGARWAL"]["scaffolds"]

ROOT = os.path.dirname(NULLS)

#write

LIBRARY_FILE = os.path.join(ROOT, f'{MER}mer.fo.pam.scaffold'+
              f".ext200.library_.tsv")
out200 = os.path.join(ROOT, f'{MER}mer.fo.pam.scaffold'+
              f".ext200_.fa")
out4096 = os.path.join(ROOT, f'{MER}mer.fo.pam.scaffold'+
              f".ext4096_.fa")

# write library file to config
config[f"{CL}.{MER}mer.{NMUTS}mut"]["200mers"] = out200
config[f"{CL}.{MER}mer.{NMUTS}mut"]["4096mers"] = out4096
config[f"{CL}.{MER}mer.{NMUTS}mut"]["library"] = LIBRARY_FILE
crw.write(config, cfn)


# In[9]:


def getNullSeqs(file, concat):
    """
    read nullomer sequence file, concatenated nullomer sequence file
    
    input
        file(str) - .txt or .tsv file to nullomer sequences
        concat(bool) - check if file is of singl enulloemrs or concatmers
        
    method
        1. make empty nulls, kmers list to collect data
        2. iterate through file lines
        3. concatmer - append null, kmer sequences
        4. single nullomers - append just null sequences. 
    return list of nullomers, matching kmers (for concatemers only)
    """
    nulls, kmers = [],[]

    # read file
    with open(file, "r") as reader:

        # iterate through lines
        for seq_id, line in enumerate(reader):
            if concat is True:
                if "nullConcat" not in line:
                    null, kmer, gcnull, gcnull, null_list = line.strip(
                    "\n").split('\t')

                    # append matched ker
                    kmers.append(kmer)
                    nulls.append(null)
            elif ">" not in line:  # the randomly sampled single nullomers. 
                null = line.strip("\n")
                nulls.append(null)
                
    return nulls, kmers

def getScaffolds(scaffold_file):

    scaffold_coor = []
    with open(scaffold_file, "r") as reader:
        for line in reader:
            coor, activity = line.strip("\n").split("\t")
            scaffold_coor.append(coor)

    print("n scaffolds", len(scaffold_coor))
    
    return scaffold_coor

def getStrand(scaffold_coor, lib, fa_dict):

    """
    return strand of the scaffold sequences

    input
        scaffold_seq_list (list)- list of scaffold coordinates
        lib (pd dataframe) - vikram's 60k common library as pandas dataframe

    method
        1. iterate through seq coordinate list
        1.1 get sequence from coordinate. 
            Because library was in hg38, but scaffolds are in hs1 (i lifted), need to harmonize on sequence using 150 bases

        2. get strand for sequence
            2.1 check for strand conflicts (where sequence matches both to a + and - strand)
        3. if strand is negative, get reverse complement (fwd), else return coor and fwd seq

    return
        scaffold_coor_fa (list) - tuple of coordinates and forward strand sequences

    """
    scaffold_coor_fa = []
    for n, coor in enumerate(scaffold_coor): 

        # get the sequence for the scaffold
        seq = fa_dict[coor]

        #2 get strand by str matching. 
        # First 150 nucleotides still return duplicate matches, 
        # but strand is consistent between matches. 

        str_ = list(set(lib.loc[lib["230nt sequence (15nt 5' adaptor - 200nt element - 15nt 3' adaptor)"].str.contains(seq[:150]), "str.hg38"]))


        # check for strand conflicts
        if len(str_)>1:
            print("discordant strands. Review sequence", seq)

        elif str_[0] == '-':
            print("\nget reverse complement", seq)
            
            fwd = "".join(reverse_complement(Seq(seq)))
            
            print("\nhere it is", fwd)
        else:
            fwd = seq
        scaffold_coor_fa.append((coor, fwd))
        
    return scaffold_coor_fa

def faToDict(fa_file):
    """
    turn all the fa coordinates into a dictionary. 
    """

    fa_dict = {}
    with open(fa_file, "r") as fasta_reader:
        for values in SimpleFastaParser(fasta_reader):
            
            name, seq = values
            fa_dict[name] = seq
            
    return fa_dict

def dissectFaId(line):
    """
    dissect mutated fa id line for kmer, nullomer information, and get the coordinate file
    """
    # process the line
    kmer, chr_, coor, strand, id_ = line.strip('\n').strip(">").split(".")
            
    # get key and value sequence to look up
    key, value_seq = kmer[:5], kmer[5:]

    bed_coor = chr_ + ":" + coor
    
    region_str = ".".join([kmer, bed_coor, strand, id_]) #TACCCATTCGAGTC.chr1.100333924-100334124.fwd
    
    
    return kmer, key, value_seq, bed_coor, region_str

def seqInsert(seq, insert_start, insert_seq):
    """
    insert sequence fragment (insert_seq) at position (insert_start) within full sequence (seq)
    
    return inserted sequence. 
    """

    insert_size = len(insert_seq)
    
    insert_end = insert_start + insert_size  # find center end

    return seq[:int(insert_start)] + insert_seq + seq[int(insert_end):]

def MPRAInsertNull(null, fa_sequence):
    
    """
    return positions where the nullomer is in fa sequence
    """
    insert_start = (len(fa_sequence)- len(null))/2
        
    null_insert = seqInsert(fa_sequence.lower(), insert_start, null.upper())  # insert the nullomer sequence into the sequence
    
    return null_insert

def mapBed(bedfile, expandedfile):

   
    bed = pd.read_csv(bedfile, sep='\t')
    bed.columns = names =["#chr", "start", "end", "strand", "id"]
    bed["coor"] = bed["#chr"] +":"+bed["start"].map(str) + "-" + bed["end"].map(str)


    ext = pd.read_csv(expandedfile, sep='\t', header=None, names =["#chr_e", "start_e", "end_e", "id"])
    ext["coor_e"] = ext["#chr_e"] +":"+ext["start_e"].map(str) + "-" + ext["end_e"].map(str)

    joint = pd.merge(ext[["id", "coor_e"]], bed[["id", "coor"]])
    
    return joint

def insert4096(null_seq, joint_bedex, fa_ex_dict, coor):

    # extend the coordinates to 4096
    extended_coor = joint_bedex.loc[joint_bedex["coor"]
                                    == coor, "coor_e"].iloc[0]

    # get the fasta file linked to the extended coordinates
    extended_seq = fa_ex_dict[extended_coor]

    # get the start to insert  # figure out where to insert the sequence
    insert_start = len(extended_seq)/2 - 100

    # insert the 200bp MPRA tile w/ nullomer variants into sequence
    null4096 = seqInsert(extended_seq.lower(), insert_start, null_seq.upper())

    return null4096, extended_seq

def prettifySeq(original, mut):
    """
    prettify sequence. 
        All matching bases will be written as "."
        All non matching bases will be written w mutated base identity. 
    """
    prettyseq = ""
    for o, m in zip(original, mut):
        if o != m:
            prettyseq += m
        else:
            prettyseq += "."

    return prettyseq

def visualInspection(kmer, null, null_ins200, inserts4096, original_seq):
    """
    print a bunch of information related to how this script is 
    inserting kmers into 200bp regions, 
    inserting those 200bp regions into 4096mers
    """
    
    print("\nkmer\n", kmer)
    
    print("\nthe nullomer sequence\n",
          null, "\n\ninserted\n", null_ins200)
    
    print("\n4096 w/ null sequence\n", inserts4096)
    
    print("\n pretty 4096", prettifySeq(
        original_seq.upper(), list(inserts4096.values())[0].upper()))

def addAdaptor(seq):
    adaptor5, adaptor3 = "AGGACCGGATCAACT", "CATTGCGTGAACCGA"
    adapted = adaptor5 + seq + adaptor3
    return adapted.upper()

def writeRow(file, row):
    """
    append line to file
    
    input
        file (str) - file to write to
        row (str) - row to write
    method
        1. open file
        2. write
        3. close
    """
    #1
    with open(file, "a") as writer:
        #2
        writer.write(row)
    #3    
    writer.close()


# # scaffold selection

# Scaffold selection? 
# 
# See: ./nullomers/bin-lock/Analysis/mpra_corr_scaffold.ipynb
#     ./nullomers/bin-lock/Sei_seqeunce_class-scaffold-minmax-selection.ipynb
# 
#     Basically, 60K sequences commonly active|inactive 
#     across 3 cell lines from Agarwal 2023 were evaluated 
#     for MPRA activity (log(rna/dna)). I quantified the mean 
#     activity across these three cell lines, then randomly sampled from  
#     the top 97.5%, IQR, and bottom 2.5% active sequences 
#     
#     Then, I predicted the sequence class from each of those 
#     common sequences. Of the sequences with interpretable predictions in Sei
#     (i.e. sei predicted active sequences as promoter/ enhancer)
#     I sampled 50 of the active and 50 of the IQR (marginal activity)
#     to use as scaffolds

# # Main

# ## set up scaffold FA resources

# In[6]:


fa_dict = faToDict(FA)
fa_ex_dict = faToDict(FA_EXT_BED)
joint_bedex = mapBed(BED, EXTENDED)
lib = pd.read_csv(JOINT_LIBRARY, header=1)


# ## get scaffolds

# In[10]:


scaffold_coor = getScaffolds(SCAFFOLDS)

# make sure all the scaffold sequences are in the forward direction
scaffold_coor_fa = getStrand(scaffold_coor, lib, fa_dict)


# In[11]:


len(scaffold_coor_fa)


# ## write single nulls

# In[ ]:


out_tsv_dict = {}  # collect data for tsv

# do single nullomers first
nulls, nullkmer = getNullSeqs(NULLS, False)
concatnulls, concatkmers = getNullSeqs(CONCAT, True)
datas = [
         ("single-null", nulls, scaffold_coor_fa),
         ("concat-null", concatnulls, scaffold_coor_fa[:1]),
         ("concat-kmer", concatkmers, scaffold_coor_fa[:1])
         ]

for datatype, seq_list, scaffold_list in datas:
    print(datatype, len(seq_list), len(scaffold_list))

    for i, seq in enumerate(seq_list):

        for j, info in enumerate(scaffold_list):
            coor, region_fa = info
            if "concat" in datatype:
                coor_ = f"concat_template-{coor}"
            else:
                coor_ = coor

        
            # get fa sequence for the coordinate, make id from coordinate and the sequence index
            region_id  = datatype + "-" + str(i)
            region_fa_id = coor_ + "_" + region_id

            # insert nulls, kmers into center of 200bp MPRA sequences.
            seq_ins200 = MPRAInsertNull(
                seq, region_fa)  # insert null concatemer

            # insert 200bp null sequence into the 4096 sequence
            seq_ins4096, original_seq = insert4096(seq_ins200, joint_bedex, fa_ex_dict, coor)

            if i == 0:  # if this is the first sequence
                
                # write the original sequences
                writeRow(out200, f">{coor}_endog_\n{region_fa}\n")
                writeRow(out4096, f">{coor}_endog_\n{original_seq}\n")

                out_tsv_dict[f"endog-{j}"] = ["endog", "endog", 
                                              i, j, coor, "+", 
                                              None, region_fa, 
                                              addAdaptor(region_fa)
                                             ]

            # write 200mers
            writeRow(out200, f">{region_fa_id}\n{seq_ins200.upper()}\n")

            # write 4096 nullomers
            writeRow(out4096, f">{region_fa_id}\n{seq_ins4096.upper()}\n")
                
            out_tsv_dict[region_fa_id] = [region_id, f"{datatype}", 
                                          i, j, coor_, "+", 
                                          seq.upper(), seq_ins200.upper(), 
                                          addAdaptor(seq_ins200)
                                         ]
            j +=1


# In[ ]:


a = [np.array(i) for i in out_tsv_dict.values()]

df = pd.DataFrame(a)
before = df.shape[0]
df= df.drop_duplicates()
print(before - df.shape[0])

df.columns=["id","exp.id", "exp.seq.id","exp.coor.id", "coor", "str", "seq_frag", "200bp_insert", "200bp_insert+adaptors"]
df["gc_frac"] = df["200bp_insert"].apply(lambda x:gc_fraction(x))
df.to_csv(LIBRARY_FILE, sep='\t', index=False)


# In[ ]:


df.shape


# # add control sequences. 

# In[ ]:


lib.head()


# In[ ]:


ctrl = lib.loc[
    (lib["category"].str.contains("negative"))|
    (lib["category"].str.contains("positive"))
     ].drop_duplicates()
ctrl = ctrl[['name',
 'category',
 'chr.hg38',
 'start.hg38',
 'stop.hg38',
 'str.hg38',
 "230nt sequence (15nt 5' adaptor - 200nt element - 15nt 3' adaptor)"]]


# In[ ]:


print(ctrl.shape)

ctrl.head()


# In[ ]:


ctrl = ctrl.fillna(-1)
ctrl["coor"] = "chr" + ctrl["chr.hg38"].map(str) + ":" + ctrl["start.hg38"].map(int).map(str) + "-" + ctrl["stop.hg38"].map(int).map(str)

ctrl.head()


# In[ ]:


df.head()


# In[ ]:


ctrl.rename(columns={"name":"id", 
                    "category":"exp.id", 
                    })


# In[ ]:


endog_seqs = list(df.loc[df["exp.id"]=="endog", "200bp_insert"])


# In[ ]:


def getStrand(seq, lib):
    str_ = lib.loc[lib["230nt sequence (15nt 5' adaptor - 200nt element - 15nt 3' adaptor)"].str.contains(seq), "str.hg38"].iloc[0]
    
    return str_


# In[ ]:


count


# In[ ]:


len(endog_seqs)


# In[ ]:




