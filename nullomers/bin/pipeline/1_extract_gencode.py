#!/usr/bin/env python
# coding: utf-8
"""
 20221122
 sarahfong
 
 1. downloads GENCODE basic gene annotations
     https://www.gencodegenes.org/human/
     
     
     This download contains the basic gene annotation on the reference chromosomes only
     This is a subset of the corresponding comprehensive annotation, including only those transcripts tagged as 'basic' in every gene
     This is the main annotation file for most users
    
2.  Use BEDOps to convert GTF -> bed

3.  Annotation filter(in config) 
    # GTF file type info: http://mblab.wustl.edu/GTF22.html
    # - CDS represents the coding sequence starting with the first translated codon and proceeding to the last translated codon. 
    # - Unlike Genbank annotation, the stop codon is not included in the CDS for the terminal exon. 
    # - The optional feature "5UTR" represents regions from the transcription start site or beginning of the known 5' UTR to the base before the start codon of the transcript. 
    # - If this region is interrupted by introns then each exon or partial exon is annotated as a separate 5UTR feature. - Similarly, "3UTR" represents regions after the stop codon and before the polyadenylation site or end of the known 3' untranslated region. 
    # - Note that the UTR features can only be used to annotate portions of mRNA genes, not non-coding RNA genes.

4.  Merge
"""


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
import zippery

config_name = os.path.join(os.getcwd(), config_tag)

config, configname = crw.read_config(config_name)


# In[6]:


# select config variables
VERSION = int(config["GENCODE"]["VERSION"])  # turn into int type
GENCODE_PATH = config["GENCODE"]["PATH"]
GENCODE_ANNOT = config["GENCODE"]["ANNOT"] # gencode annotation to filter on 


# # functions 

# In[7]:


def download_gencode(path, version):
    """
    
    Returns local gencode files 
        if not local - downloads metagene, basic.annotation.gtf files from gencode FTP
    
    inputs 
        path (str) - name to local gencode datapath
        version (int) - gencode version to download
        
        
    method
        1. make FTP version path for gencode
        2. make local gencode file names
        3. check if basic annotation version file exists locally
            3a - if not, download 
        4. check if metagene version file exists locally
            4a - if not, download
        
        return local copies of the download files
    """
    
    #1 FTP path
    FTP_PATH = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{version}"    
    
    
    #2 Local file names
    basic_gene_file=os.path.join(path, f"gencode.v{version}.basic.annotation.gtf.gz")
    meta_gene_file=os.path.join(path, f"gencode.v{version}.metadata.HGNC.gz")

    #3 BASIC GENE ANNOTATION .GTF file
    if os.path.exists(basic_gene_file) is False and os.path.exists(basic_gene_file.strip(".gz")) is False:
        cmd = f"wget {FTP_PATH}/gencode.v{version}.basic.annotation.gtf.gz"
        
        print("downloading basic gene annotations (ref chromosome only)", version)
        os.system(cmd)
    else:
        print("already downloaded basic gene annotation (ref chromosome only)", version)
        
    #4 METAGENE FILE - gene names
    if os.path.exists(meta_gene_file) is False:
        cmd = f"wget {FTP_PATH}/gencode.v{version}.metadata.HGNC.gz"
        
        print("downloading meta data", version)
        os.system(cmd)
    else:
        print("already downloaded metadata", version)
        
    return basic_gene_file, meta_gene_file


# In[8]:


def gtf2bed(gtf_file):
    """
    convert gtf2bed using bedops tool
    
    input
        gtf_file (str) with full path to the .gtf file
    
    method
        1. create an output .bed file
        2. compile the commandline command. 
        3. run through the command line
    
    return
        output .bed file (str)
    """
    
    #1
    out_bed = gtf_file.strip(".gtf") + ".bed"
    
    #2
    cmd = f"gtf2bed < {gtf_file} > {out_bed}"
    
    #3
    if os.path.exists(out_bed) is False:

        print("converting to bed\n\n", cmd)
        subprocess.call(cmd, shell=True)
        
    elif os.path.getsize(out_bed)==0:
        print("converting to bed\n\n", cmd)
        subprocess.call(cmd, shell=True)

    else:
        print("already made .bed")
    
    return out_bed


# In[9]:


def filterGencodeBed(bed_file, annotation):
    """
    filters gencodebed.bed file on annotations
    removes chrM annotations
    
    input 
        bed_file (str) - full path to the gencode bed file
        annotation (str) - annotation value to filter on
        
    method
        1. open the dataframe in pandas
        2. remove chrM annotations
        3. keep only annotation matches



    annotation possibilities:
        'exon', 
        'gene', 
        'transcript', 
        'UTR', 
        'start_codon', 
        'CDS', 
        'stop_codon', 
        'Selenocysteine'
        
        CDS = first translated codon + last translated codon. No terminal codon.
        Exon = sequence that is included in mature mRNA, may or may not get translated
        5'UTR = TSS to transcript-1 position
        3'UTR = stop codon+1 to poly-a, or end of 3'UTR transcript
    """
    
    #1
    df = pd.read_csv(bed_file, sep ='\t', header=None)
    
    beforesize = df.shape[0]
    
    #2 remove chrM
    df = df.loc[df[0]!= "chrM"]
    
    aftersize = df.shape[0]

    print("removing chrM annotations, n=", beforesize-aftersize)

    #3 keep only gene annotations
    filtered_df = df.loc[df[7]==annotation]
    print("keeping only annotations", annotation, filtered_df.shape)
    
    
    return filtered_df


def flattenGencodeBed(bed_file):
    
    """
    Return a GENCODE bed file with NO OVERLAPPING coordinates. 

    Goal - want to shuffle into background where coordinates are unique. 
        This will reduce bias from genes with many exons identified from different sources. 
        Not flattening, inflates representation of exons w/ same or similar coordinates annotated more than once in GENCODE.
    
    input 
        bed_file (str) - full path to bedfile
    
    method
        1. make FLAT_GENCODE.bed file name
        2. if .bed does not exist, 
            2a. use pybedtools merge function (any coodinate overlapping 1 bp will be merged)
    
    return
        FLAT_GENCODE (str) - path to flattend file. 
            
    """
    FLAT_GENCODE = bed_file.strip(".bed") + "-merged.bed"
    
    if os.path.exists(FLAT_GENCODE) is False:
        gencode = pbt.BedTool(bed_file)

        # merge all elements that overlap 1bp
        flattend = gencode.merge(output=FLAT_GENCODE)
    
    else:
        print("you flattened this already")
    
    return FLAT_GENCODE


def main(argv):


    # ## download gencode annotations

    # files to write
    BED = os.path.join(GENCODE_PATH, f"gencode.v{VERSION}.basic.annotation.bed")
    ANNOT_BED = os.path.join(GENCODE_PATH, f"gencode.v{VERSION}.basic.annotation-{GENCODE_ANNOT}.bed")

    if os.path.exists(BED) is False:

        basic_gtf, meta_gtf = download_gencode(GENCODE_PATH, VERSION)

        # unzip the gtf file
        unzipped_gtf = zippery.unzip_file(basic_gtf) 


    # ## gtf2bed 
    # https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gtf2bed.html

    # In[ ]:


    # convert gtf 2 bed w/ BEDOps - 

    if os.path.exists(BED) is False:

        gtf_bed = gtf2bed(unzipped_gtf)


    # rezip the gtf file

    if os.path.exists(BED) is False:

        zippery.rezip_file(unzipped_gtf) 


    # ## filter for annotation


    filtered_annot = filterGencodeBed(BED, GENCODE_ANNOT)


    # write
    
    filtered_annot.to_csv(ANNOT_BED, sep='\t', header=None, index=False)
    
    # bedtools merge 
    FLAT_GENCODE = flattenGencodeBed(ANNOT_BED)

    # write files to config
    
    config["GENCODE"]["BED"] = BED
    config["GENCODE"][f"{GENCODE_ANNOT}_BED"] = ANNOT_BED
    config["GENCODE"]["MERGED"] = FLAT_GENCODE

    crw.write_config(config, configname)
    
if __name__ == "__main__":
    main(sys.argv[1:])
