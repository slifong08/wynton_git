import glob
import os
import pybedtools as pbt
import subprocess
import sys
 
import numpy as np
import pandas as pd


# append path
sys.path.append("/wynton/home/ahituv/fongsl/tools/py_/")

# import config reader
import config_readwrite as crw
import zippery

config_tag, NBASES, ARRAY_IND = sys.argv[1], sys.argv[2], sys.argv[3]

NBASES = int(NBASES)

#config_tag, NBASES, ARRAY_IND = "config-exon.ini", "500", "3"
print(ARRAY_IND, type(ARRAY_IND), config_tag, NBASES, type(NBASES))

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

# get QSUB array
ARRAY = config["QSUB"]["ARRAY"]  # tab separated file. 

# linsight files, wiggle track executables
LINSIGHT_BW = config["LINSIGHT"]["LINSIGHT_BW_HG38"]

LINSIGHT_PATH = config["LINSIGHT"]["PATH"]
SRC_BW = config["SRC"]["bigwigsummary"]


config["LINSIGHT"]["FLANK_BP"]= str(NBASES)


# # functions
def getArrayMatch(arrayfile, array_ind):

    RETURN_BED = ""
    ### return .bed file from array file. 
    with open(arrayfile, "r") as array:
        for line in array.readlines():

            JOB_NUM, BEDFILE = (line.strip("\n")).split("\t")
            if str(JOB_NUM) == str(array_ind):  # get the number corresponding to the file. 
                print("match array",JOB_NUM, BEDFILE)
                RETURN_BED = BEDFILE
                
                
    return RETURN_BED
    

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

    #4
    with open(bedfile, "r") as file:
        for line in file:
            print("original", line)
            break  # only print the first original line in the file
            
    with open(out, "r") as ex:
        for line in ex:
            print("expanded", line) 
            break  # only print he first expanded line in the file
    
    return out

# ## parse bigwig for bed coordinates


def getBwBedLineSummary(bw, bed_line_list, nbases, outfile, src):
    
    """
    write bed coordinates + linsight values for a bed-coordinate region. 
    
    input
        bw (str) - path to bigwig file to parse
        bed_line_list (list) - list with CHR, START, STOP, and 4th column annotations. 
            Should correspond to 1 .bed line.
        nsteps (int) - steps to summarize the linsight sequence in. 
            e.g. If 1000 bases in sequences and nsteps == 1000, will return each bp linsight score. 
            e.g. if 1000 bases in sequences and nsteps == 10, will average (?) 100bp linsight scores into 1 output. 
        outfile (str) - full path to write file. 
        src (str) - path to ./bigWigSummary executable. 
        
    method
        1. unpack bed_line_list into CHR, START, STOP, and sample ID variables
        2. run bigWigSummary command, get phyloP scores for the bed coordinates (output). 
        3. make new line (str) w/ bedcoordinates, id, linsight results
        4. write bed coordinates and linsight to the outfile
        
        
    """
    
    #1
    chr_, start_, end_, nid = bed_line_list
    nsteps = int(nbases)*2
    
    # TO DO: insert a check for existing file, number of lines in file. Concerned about parallelized runs and writing to wrong coordinates. 
    
    #2              
    cmd = f"{src} {bw} {chr_} {str(start_)} {str(end_)} {str(nsteps)}"
    #print(cmd, "\n")
 
    result = subprocess.run([src, bw, chr_, str(start_), str(end_), str(nsteps)], stdout=subprocess.PIPE)
 
    #3 make a str (tab-delimited) of the bed coordinates and linsight results
    newline = "\t".join(bed_line_list) +"\t"+ result.stdout.decode('utf-8') +"\n"
    #print(newline)
    
    #4 
    with open(outfile, "a") as writer:
        writer.write(newline)
        writer.close()

# get summary stats for linsight scores

# get summary stats for phylop scores

def q025(x):
    return x.quantile(0.025)

def q975(x):
    return x.quantile(0.975)


def getSummaryStats(file, value_col):
    """
    return summary stats for each position in a bed element
    
    input
        file (str) - path to extracted bw vector file from above run
        value_col (str) - name of column to calculate values on. 
        
    method
        1. make a pandas dataframe out of the dictionary
            
            col 0-3 = CHR, START, STOP, Nullomer ID
            col 4-1003 = each bp position score flanking nullomer
            
            col 503 = pre-nullomer locus.
            
            1.1 Formatting - SKIP BAD LINES
                - drop_duplicates() 
                - dropna()

    
        2. turn dataframe into long form data, where 1 column is the position, and another column is the value
            index on CHR, START, STOP, Nullomer ID
            
            - DROP NAN values. 
            
        3. per position, bootstrap confidence intervals. 
        
        4. compute summary stats per position across nullomer sequence for linsight scores.
            4.1 also correct for
        
        5. write summary stats to outfile
    
    return 
    
        summary stats outfile (str) - write summary stats as file, 
            ## so that stats can be done on empirical summary statistics. 
            
        
    """
    out = file.strip(".txt") + "-SUMMARY_STATS.txt"
    zipped = file + ".gz" 

    if os.path.exists(out) is False:
        
        # if re-running stats
        if os.path.exists(zipped) is True:
            zippery.unzip_file(zipped)
            
        
        #1 make a pandas dataframe
        df = pd.read_csv(file, 
                         sep = '\t', 
                         header = None, 
                         low_memory = False, 
                         error_bad_lines=False).drop_duplicates().dropna() #1.1

        #2

        melted = pd.melt(df, id_vars= df.columns[:4], var_name="pos", value_name=value_col)
        print(melted.shape, melted.head())

        melted[value_col] = melted[value_col].astype(float)  # change data type
        

        f = {value_col: [q025, q975]}
        ci = melted.groupby('pos').agg(f).reset_index()

        #4
        summary_stats = melted.groupby("pos")[value_col].describe().reset_index()
        summary_stats= summary_stats.join(ci)

        print(list(summary_stats))
        #summary_stats.drop( "('pos', '')", axis=1)

        summary_stats["pos"] = summary_stats["pos"] - 4 # correct bp position 
                            

        #5 
        summary_stats.to_csv(out, sep ='\t', index=False)

        zippery.rezip_file(file)
    
    return out

# # Main 

# ## exonic, non-exonic; expand; extract bw values

def main(argv):
    
    # combined function that expands mutations and extracts linsight scores
    #out = expandAndBwExtract(BED, NBASES, LINSIGHT_BW, LINSIGHT_PATH)
    """
    extract and summarize linsight for expanded bed coordinates. 
    
    inputs
        BED (str) - bed file name w/ full path from ARRAY
        LINSIGHT_PATH (str) - full path to output file directory
        NBASES (int) - number of bases to expand on each side
        LINSIGHT_BW (str) - full path to linsight bigWig file. 
        
        
    method
    0. makes an outfile. 
        0.1 if scores have not been extracted already
        0.2 if scores have been and the file is zipped
    1. takes locus and expands the bed coordinates of the locus. 
    2. extracts the linsight base scores for the expanded locus, writing the scores as a vector
    3. write summary stats of file. (Summary stats at each bp)


    return 
        outfile (str) - name to output values. 
    """
    
    BED = getArrayMatch(ARRAY, ARRAY_IND)

    #0 make the outfile
    name = (BED.split("/")[-1]).strip(".bed") + "-linsight.txt"
    outfile = os.path.join(LINSIGHT_PATH, name)
    print("MAKE THIS FILE FIRST", outfile)
         
    
    if os.path.exists(outfile +".gz") is False and os.path.exists(outfile) is False:

        #1 expand bed coordinates
        print("expanding", BED, "n bases", NBASES)

        outex = expand_bedcoor(BED, int(NBASES))

        print("EXPANDED", outex)

        #2 get linsight values from bigwig and write to the outfile. 
        with open(outex, "r") as ex:

            print("reading expanded file", outex)

            for line in ex.readlines():
                bedlist = (line.strip("\n")).split("\t")
                getBwBedLineSummary(LINSIGHT_BW, bedlist, NBASES, outfile, SRC_BW)

        print("wrote linsight scores to", outfile)

    # 3
    summary_outfile = getSummaryStats(outfile, "linsight")

    print(summary_outfile)
    
    # # write config

    crw.write_config(config, configname)

if __name__ == "__main__":
    main(sys.argv[1:])