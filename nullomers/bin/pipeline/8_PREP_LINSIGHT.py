import glob
import os
import subprocess
import sys

# append path
sys.path.append("/wynton/home/ahituv/fongsl/tools/py_/")

# import config reader
import config_readwrite as crw
import zippery

# INPUTS
config_tag = sys.argv[1]

# CONFIG
config, configfile_name = crw.read_config(os.path.join(os.getcwd(), config_tag))

DATAPATH = config["DATA"]["PATH"]

# LINSIGHT resources
LINSIGHT_BW_URL = "http://compgen.cshl.edu/LINSIGHT/LINSIGHT.bw"  # precomputed scores. In hg19. 
LINSIGHT_PATH = os.path.join(DATAPATH, "linsight")
LINSIGHT_BW = os.path.join(LINSIGHT_PATH, "LINSIGHT.bw")
LINSIGHT_BED = os.path.join(LINSIGHT_PATH, "LINSIGHT.hg19.bed")
LINSIGHT_BED_HG38 = os.path.join(LINSIGHT_PATH, "LINSIGHT.liftOver.to.Hg38.bed")
LINSIGHT_BW_HG38 = os.path.join(LINSIGHT_PATH, "LINSIGHT.liftOver.to.Hg38.bw")

# bigwig to bedgraph executable
BIGWIGTOBEDGRAPH_EXE = "/wynton/home/ahituv/fongsl/nullomers/src/bigWigToBedGraph"
BEDGRAPHTOBIGWIG_EXE = "/wynton/home/ahituv/fongsl/nullomers/src/bedGraphToBigWig"
CHR_SIZE = "/wynton/home/ahituv/fongsl/dna/hg38/hg38.chrom.sizes"

# LIFTOVER 
# local liftover executable + files. 
LIFTOVER_EXE = "/wynton/home/ahituv/fongsl/nullomers/src/liftOver"
LIFTOVER_PY = "/wynton/home/ahituv/fongsl/tools/evo/liftover_bed-wynton.py"
LOCAL_CHAIN_PATH = "/wynton/home/ahituv/fongsl/dna"


FROM = "hg19"
TO = "Hg38"
CHAIN_URL = f"https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/{FROM}To{TO}.over.chain.gz"
LOCAL_CHAIN = os.path.join(LOCAL_CHAIN_PATH,f'{FROM}To{TO}.over.chain.gz')

###
# LINSIGHT configuration
###

def configureLinsight(linsight_path, linsight_bw, linsight_bw_url, bigwigtobedgraph_exe, linsight_bed, config):
    
    """
    configure linsight .bed file
        - need to convert to .bed in order to liftOver
        
    inputs
        linsight_path (str) - path to linsight data (originally in hg19)
        linsight_bw (str) - linsight.bw file, full path, originally downloaded from LINSIGHT_BW_URL
        linsight_bw_url (str) - web path to original linsight.bw file from the Siepel lab.
        bigwigtobedgraph_exe (str) - local path to bigWigToBedGraph exectuable from UCSC
        linsight_bed (str) - file to write linsight.bed data. 
        config (configparser object) - config file to write inputs/outputs to
        
    method
        1. check to make sure datapath exists, if not, make directory
        2. add linsight section to config
        3. download local linsight.bw, if not already downloaded. 
        4. .bw --> .bed for all hg19 linsight pre-computed estimates. 
        
    return 
        none. 
        
    """
    #1 if the datafolder does not exist, make it

    if os.path.exists(linsight_path) is False:
        os.mkdir(linsight_path)
        print("made LINSIGHT data folder")


    #2 make a section in config file for LINSIGHT. 
    crw.check_section(config, "LINSIGHT")

    #3 download linsight file
    zipped_bw = linsight_bw + ".gz"
    
    if os.path.exists(linsight_bw) is False and os.path.exists(zipped_bw) is False:
        print("downloading LINSIGHT.bw (hg19!) from internet")
        os.system(f"wget -O {linsight_bw} {linsight_bw_url}")

    #4 turn linsight BW -> BED
    zipped_bed = linsight_bed + ".gz" # output file
    
    if os.path.exists(linsight_bed) is False and os.path.exists(zipped_bed) is False:
        
        print("converting BW -> BED")
        os.system(f"{bigwigtobedgraph_exe} {linsight_bw} {linsight_bed}")
        zippery.rezip_file(linsight_bw) # rezip the linsight bw file
        
        
    
###
# Liftover configuration
###

def configureLiftOver(local_chain, chain_url, from_, to_, liftover_py, linsight_bed_to, linsight_bed_from, config):
    
    """
    configure liftOver chains, liftOver linsight.bed file to hg38
    
    inputs
        local_chain (str) - path to local chain file for liftover
        chain_url (str) - web address to chain url from golden path. 
        from (str) - genome build to lift from 
            *** NOTE - ALL LOWER ***
        to (str) - genome build to lift to 
            *** NOTE - FIRST LETTER MUST BE UPPER ***
        liftover_py (str) - path to sarah's custom script that liftOvers 
        linsight_bed_to (str) - path to linsight_bed in hg19 to be liftedOver
        linsight_bed_from (str) - path to file that will be written in liftOver
        config (configparser obj) - config to write inputs and outputs to. 
        
    method
        1. make a section in the config file for liftover inputs and outputs
        2. download chain files from golden path if not already downloaded locally 
        3. liftover entire linsight file from hg19 --> hg38
        
    return 
        None
    """

    #1 make a section in config file for liftover. 
    crw.check_section(config, "LIFTOVER")
    
    #2 download chain file. 
    if os.path.exists(local_chain) is False:

        print(f"downloading chain file(from {from_}, to {to_}) from internet")
        os.system(f"wget -O {local_chain} {chain_url}")

    #3 do liftover hg19 --> hg38
    run_liftover = f"python {liftover_py} -b {linsight_bed_from} -f {from_} -t {to_}"
    
    zipped_tobed = linsight_bed_to + ".gz"
    zipped_frombed= linsight_bed_from + ".gz"
    
    if os.path.exists(linsight_bed_to) is False and os.path.exists(zipped_tobed) is False:
        
        if os.path.exists(zipped_frombed) is True: # check if from bed is unzipped
            zippery.unzip_file(zipped_frombed)
            
        print(run_liftover)
        os.system(run_liftover)
        
        zippery.rezip_file(linsight_bed_from)# rezip the file
        
        # sort liftover bedfile
        sortBed(linsight_bed_to)

    elif os.path.exists(zipped_tobed) is True:
        zippery.unzip_file(zipped_tobed) # unzip the file
        
        print(run_liftover)
        os.system(run_liftover)
        
        zippery.rezip_file(linsight_bed_from)# rezip the file

        # sort liftover bedfile
        sortBed(linsight_bed_to)
        

def bedToBigWig(linsight_bed, linsight_lifted_bw, bedgraphtobigwig_exe, chr_size):
    
    """
    .bed --> .bw
    
    input 
        linsight_bed (str) - full path to liftedover linsight bed file. FOUR columns only. 
        linsight_lifted_bw (str) - full path to liftedover linsight bw file (to write)
        bedgraphtobigwig_exe (str) - full path to bedgraphtobigwig executable from ucsc.
        chr_size (str) - full path to lifted over genome's chromosome sizes
        
    method
        1. merge liftover file. Take MEDIAN of LINSIGHT scores. Merged elements must overlap by at least 1 bp (-d -1)
        2. if linsight_lifted_bw does not exist, run bedToBigWig
        
    return
        None
        
    """
    #1
    merged = linsight_bed.strip(".bed") + ".merge.bed"

    if os.path.exists(merged) is False:
        print("merging lifted")

        # bedtools merge
        os.system(f"bedtools merge -i {linsight_bed} -d -1 -c 4 -o median > {merged}")

        zippery.rezip_file(linsight_bed) # zip the old file
        
    if os.path.exists(linsight_lifted_bw) is False:
        
        # bedGraphToBigWig command
        cmd = f"{bedgraphtobigwig_exe} {merged} {chr_size} {linsight_lifted_bw}"
        print(cmd)
        os.system(cmd)
    else:
        print("converted bed back to bw already")
        
    return merged
        
        
def addLinsightConfig(config):
    """
    write inputs and outputs to config LINSIGHT section
    """
    
    # add to linsight config
    config["LINSIGHT"]["LINSIGHT_BW_URL"]= LINSIGHT_BW_URL
    config["LINSIGHT"]["PATH"]= LINSIGHT_PATH
    config["LINSIGHT"]["LINSIGHT_BW"]= LINSIGHT_BW
    config["LINSIGHT"]["LINSIGHT_BED_HG19"]= LINSIGHT_BED
    config["LINSIGHT"]["LINSIGHT_BED_HG38"]= LINSIGHT_BED_HG38
    config["LINSIGHT"]["LINSIGHT_BW_HG38"] = LINSIGHT_BW_HG38
    

def addLiftOverConfig(config):
    """
    write inputs and outputs to config LIFTOVER, SRC sections
    """
    
    # add to liftover config

    config["LIFTOVER"]["LIFTOVER_PY"]=LIFTOVER_PY
    config["LIFTOVER"]["LOCAL_CHAIN_PATH"]=LOCAL_CHAIN_PATH
    config["LIFTOVER"]["FROM"]=FROM
    config["LIFTOVER"]["TO"]=TO
    config["LIFTOVER"]["CHAIN_URL"]=CHAIN_URL
    config["LIFTOVER"]["LOCAL_CHAIN"]=LOCAL_CHAIN

    # add to src config
    config["SRC"]["LIFTOVER_EXE"]=LIFTOVER_EXE
    config["SRC"]["BIGWIGTOBEDGRAPH_EXE"]=BIGWIGTOBEDGRAPH_EXE
    config["SRC"]["BEDGRAPHTOBIGWIG_EXE"]=BEDGRAPHTOBIGWIG_EXE
    config["SRC"]["CHR_SIZE"]=CHR_SIZE



# intersection of non-coding w/ linsight

def intersectBedxLINSIGHT(bed, linsight_path, linsight_bed):
    
    """
    intersect bed file w/ linsight file
    
    inputs
        bed (str) - bed file to instersect w linsight file
        linsight_path (str) - full path to directory to write file
        linsight_bed (str) - linsight bed file w/ full path
 
    
    method
        1. make a result file handle w. full path using bed file as input. 
        2. if the result file does not exist, do linsight intersection
            -wa : write all bed file rows that overlap linsight
            -wb : write corresponding linsight values for each bed. 
        
    return
        resulting bed 
    """
    
    filename = bed.split("/")[-1].strip(".bed") + "xLINSIGHT.bed"
    out = os.path.join(linsight_path, filename)
    
    if os.path.exists(out) is False:
        cmd = f'bedtools intersect -a {bed} -b {linsight_bed} -wa -wb > {out}'
        print("intersecting bed w/ linsight", out)
        os.system(cmd)
    else:
        print("intersected")
        
    
    return out


def sortBed(file):
    os.system(f"bedtools sort -i {file} > t && mv t {file}")

    
## Main

def main(argv):
    
    """
    return intersection of LINSIGHT scores
    """
    
    # configure LINSIGHT
    configureLinsight(LINSIGHT_PATH, LINSIGHT_BW, LINSIGHT_BW_URL, BIGWIGTOBEDGRAPH_EXE, LINSIGHT_BED, config)
    
    # configure LiftOver for entire LINSIGHT dataset
    configureLiftOver(LOCAL_CHAIN, CHAIN_URL, FROM, TO, LIFTOVER_PY, LINSIGHT_BED_HG38, LINSIGHT_BED, config)
    
    
    # convert lifted .bed -> .bw
    merged = bedToBigWig(LINSIGHT_BED_HG38, LINSIGHT_BW_HG38, BEDGRAPHTOBIGWIG_EXE, CHR_SIZE)
    
    
    # add LINSIGHT files to config
    addLinsightConfig(config)
    
    # add LIFTOVER and SRC files to config
    addLiftOverConfig(config)
    
    config["LIFTOVER"]["LINSIGHT_BED_HG38"]=merged  # replace the regular bed w/ the merged bed. 

    # write to config 
    crw.write_config(config, configfile_name)
    

    # do linsight intersection in .bed format
    #FILExLINSIGHT = intersectBedxLINSIGHT(BED_EXP, LINSIGHT_PATH, LINSIGHT_BED_HG38)

if __name__ == "__main__":
    main(sys.argv[1:])
    