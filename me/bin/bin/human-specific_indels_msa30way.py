#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import AlignIO
import os, sys
import numpy as np

INDEX = sys.argv[1]

# append path""
sys.path.append("/wynton/home/ahituv/fongsl/tools/py_/")

# import config reader
import config_readwrite as crw
import zippery
from chr_functions import make_chr_list


# source path on wynton. Don't manipulate files here. 
#INDEX = 19


# # functions 



def copySrcToLocalMaf(index):
    """
    make local MAF copy for index
    #  Get source files.  
    
    input 
        index (str) - index number corresponding to chr_list index, which has CHRi identity
        
    method
        1. make list of chr autosomes by calling make_chr_list function
        2. append sex chromosomes to list. 
        3. retreive chromosome number and zipped maf file handle
        4. instantiate variables w full paths to src and local data
        5. copy file to local path (if not already there, already unzipped)

    return
        CHR (str) - identity of chromosome
        DATA_PATH (str) - path to output data
    
    """
    
    #1 make the chromosome list. 
    chr_list = make_chr_list()
    
    #2
    chr_list.extend(['chrX', 'chrY'])
    
    #3 

    CHR = chr_list[int(index)-1] # minus 1 to correct to zero-base (qsub is 1 based) 
    TEST= f"{CHR}.maf.gz"

    
    #4
    MAF_PATH="/wynton/group/databases/goldenPath/hg38/multiz30way/maf" #source path on wynton. Don't manipulate files here. 
    DATA_PATH = "/wynton/home/ahituv/fongsl/dna/hg38/multiz30way"

    MAF_SRC = os.path.join(MAF_PATH, TEST)
    MAF_LOCAL = os.path.join(DATA_PATH, TEST)

    #5 ## move source files to local directory 

    os.chdir(DATA_PATH)  # change directory to local path

    #6 instead, move files to local path
    if os.path.exists(MAF_LOCAL) is False and os.path.exists(MAF_LOCAL.strip(".gz")) is False:
        os.system(f"cp {MAF_SRC} {DATA_PATH}")


    return CHR, DATA_PATH

def unzipLocalMaf(chr_, data_path):
    """
    unzip local maf chromosome number corresponding to index

    
    input 
        chr_ (str) - chr_ to unzip
        data_path (str) - path to dir
        
    method
        1. make str to chr maf
        2. change directories to data_path (did this once, but do it again
        3. if local maf is zipped, unzip

    return
        none
    
    """
    
    #1
    TEST= f"{chr_}.maf.gz"
    print(TEST)
   
    MAF_LOCAL = os.path.join(data_path, TEST)

    #2 
    os.chdir(data_path)  # change directory to local path

    #3 unzip
    if os.path.exists(MAF_LOCAL.strip(".gz")) is False:
        zippery.unzip_file(TEST)  # unzip file


def parseDelsIns(align_array, species_vector, primate_block_index, maf_block_start):
    
    """
    return a list of human-gaps and human-insertions
        human-gaps = gap in humans, but no other species at one or multiple alignment block positions
        human-insertions = insertion in humans, but all other species have gap at one or multiple alignment block positions
    
    input
        align_array (np.array) - array of multiple sequence alignment block object from AlignIO
        species_vector (list) - list of species ordered in maf block
        primate_block_index(list) - indexes of the closest related primate genomes to consider. 
    
    method
        1. instantiate gap, ins lists to collect gap, insertion positions
        2. parse through each locus
        3. get bases at locus
            3.1 get the human base (str), change CHAR to upper (lower = repeats in alignment)
            3.2 bases from other species (np.array) at that single locus, 
                keep rows corresponding to primate block, 
                change CHAR to upper (lower = repeats in alignment)
            
        4. conservation - check whether the human base is found in any of the other species at this locus. 
        5. gaps - check whether the human base is a gap AND if it is a gap in any other species at this locus. 
            5.1. append human gaps to gap list
        6. insertions - check whether the human base is a nucleotide AND if all other species have a gap. 
            6.1. append human insertion to the insertion list
            
        7. substitutions - human variant cannot be the same as any of the other primates. 
            
    return
        gap (list) - list of positions with human-specific gap (and gap in no other species. False negative issue if gap is spuriously assigned beacuse of poor sequence covergae in other species)
        ins (list) - list of positions with human-specific insertion (and only gap in other species)
    """
    
    # 1 lists to collect the positions w/ gaps, insertions specific to humans
    current_genome_pos = maf_block_start - 1
    gap, ins, sub = [], [], []  
    
    # 2 parse through each base
    for pos in np.arange(align_array.shape[1]):
        current_genome_pos+=1
        
        #print(pos, current_genome_pos)

        #3.1 get the human base
        hu = align_array[0, pos]
        
        #3.2 get the other species(') base(s) at primate indexes
        other_vals = align_array[:, pos][primate_block_index]
        
        #print(set(hu), set(other_vals))

        # 4 check for conservation of the human base
        if hu == "-":
            current_genome_pos-=1 # correct for position in the human genome
            
        if hu in other_vals:  # MOST COMMON CASE. If the human base is found in the other species vector, 
            continue  # the human base has been observed in another species, so we move on

        # 5 test for human gaps
        elif hu == "-" and "-" not in set(other_vals): # gap in human (i.e. human base == "-", and there is no gap in other species vector. )

            gap.append(current_genome_pos)
            
            #print( "\nhugap - pos", pos, align_array[:, pos], align_array.shape, species_vector)

        # 6 test for human insertions
        elif len(set(other_vals)) == 1 and "-" in set(other_vals):  # this base is observed in humans, but other species have a gap AND ONLY A GAP
            
            ins.append(current_genome_pos)
            
            #print( "\nhu ins - pos", pos, align_array[:, pos], align_array.shape, species_vector)
        
        #7 test for human-specific variants
        elif set(hu).difference(set(other_vals))==set(hu) and len(set(other_vals))==1:

            sub.append(current_genome_pos)
            
            #print( "\nhu sub - pos", pos, align_array[:, pos], align_array.shape, species_vector)
            
    return gap, ins, sub





def getConsecCoorIntervals(input_list_, mutation_type):
    
    """
    Return a list of gap/insertion genome coordinates
    
    input
        input_list_ (list) - list of zero-based coordinates where gaps or deletions live
        mutation_type (str) - INS or DEL

    method
        1. copy input list and add a "-1" to the copied list. 
            "-1" will ensure that every value in list gets processed by rules
            
        2. deletion parsing. In hg38, human gaps do not have a genomic locus address. 
            Therefore, the coordinates will be start and start +1 where human genome aligns.
            I count the number of gap bases in a region and record this value. Should match genome browser and n deleted bases.
            
        3. insertion parsing. Count the number of consecutive hg38 bases (nbases) that have no alignment in the other closest primate-species. 
        
        4. Substitution parsing. Write each sub as zero-base .bed coordinates. Each sub should be 1 bp long.
            
        
    return
        genomic_coord (list) - list of THREEPUL (genomic_Start, genomic_stop, nbases in gap/indel) for each interval in the input_list
        
    """
    
    #1
    coord = []
    
    # copy the input list
    input_list = input_list_.copy()
    
    # append a -1 to stop the run 
    input_list.append(-1)
    
    #2 deal with deletions
    if mutation_type == "DEL":
        
        # set vars for counting genomic position, n gaps
        last_start, nbases = 0, 1
        
        # parse through each listed position
        for i, start in enumerate(input_list):
            
            # handle the first index
            if i==0:
                # reassign last start to start position
                last_start = start
                #print(start, "start!")
                
            # if the first index is the last in
            elif input_list[i] == input_list[i-1]:
                nbases +=1
                #print(nbases, "running a stretch")
                
            # if new start, or last index in list
            elif input_list[i] != input_list[i-1]:
                
                # settle the last stretch of coordinates
                last_coordinates = (last_start, last_start+1, nbases)
                #print(last_coordinates, "last coor")
                
                #append coordinates
                coord.append(last_coordinates)
                
                if i != input_list[-1]:
                    #reset counter 
                    last_start, nbases = start, 1
                    #print(last_start, start, nbases, "new start")
    
    # deal with insertions
    elif mutation_type == "INS":
        last_start, nbases = 0, 1
        
        for i, start in enumerate(input_list):
            
            # deal with start of list
            if i ==0:
                # reassign last start to start position
                last_start = start
            
            # if start and last start are 1bp apart
            elif start == (last_start+1):
                nbases +=1
                
            # if new start, or last index in list
            elif input_list[i] != (input_list[i-1] +1):
                
                # settle the last stretch of coordinates
                last_coordinates = (last_start, last_start+nbases, nbases)
                #  append last coordinates
                coord.append(last_coordinates)
                
                if i != len(input_list):
                    
                    #reset counter 
                    last_start, nbases = start, 1
                    #print(last_start, start, nbases, "new start")
                    
    # deal w subs                        
    elif mutation_type == "SUB":
        for i, start in enumerate(input_list):
            if input_list[i] != -1:
                coordinates = (start, (start+1), 1) 
                coord.append(coordinates)
            

    return coord



def writeInfo(chr_, species_vector, nspecies, maf_block_len, 
              maf_block_start, genomic_indel_list, mutation_type, outfile):
    
    """
    write .bed file with genomic coordinates of indels, number of species in maf block alignment, maf_block_len
    
    inputs
        chr_ (str) - chromosome
        species_vector (list) - list of species in maf block
        nspecies (int) - number of species in maf block
        maf_block_len (int) - length of the maf block
        maf_block_start (int) - hg38 start of maf block
        genomic_indel_list (list) - list of tuples w/ start, stop coordinates of gaps or insertions. 
        outfile (str) - file to write to, with full path please!
        
    method
        1. instantiate a lines list to collect lines
        2. write line with info for each genomic interval and maf block
        3. open the outfile to write IN APPEND MODE. DO NOT RUN MULTIPLE TIMES WITHOUT CHECKING FIRST
        4. write the lines. 
        5. close the file. 

    return
        nothing to return. 
    """
    
    lines = []  # make a line for each genomic interval in genomic_indel_list
    
    for start, stop, nbases in genomic_indel_list:
        
        species_str = ",".join(species_vector) # turn vector into a strin
        
        line = "\t".join([
            chr_, str(start), str(stop), 
            str(nspecies), str(maf_block_start), str(maf_block_len), 
            species_str, 
            mutation_type, str(nbases), 
        ]) + "\n"
        lines.append(line)
        
    with open(outfile, "a") as writer:
        
        for line in lines:
            writer.write(line)
        
        writer.close()      


# # Main

# # Deletion / insertion definition
# - hg38
# 
# if hg38 gap/insertion is not in >=2 of the following species:
# 
# - gorGor5 (2016 SMRT), 
# - rheMac8 (2015), 
# - panTro5 (2016)
# - panPan2 (2015)
# - DNU - ponAbe2 (reciprocal best, not syntenic, from 2007)
# - Other possible: Gibbon (2012) (GGSC Nleu3.0/nomLeu3)
# 
# - https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz30way/



def main(argv):


    """
    Write human-specific substitutions, insertions, and deletions per input chromosome. 
    
    
    Input 
        chr (str) - chromosome number
        data_path (str) - full path to folder to dump data. 

    Method
        PREPROCESSING - move .maf.gz from wynton source to local dir. 

        0. instantiate outfiles to write. one for insertions, one for gaps
        1. unzip and parse the maf file
        2. turn each blcok into a numpy array, change all lower case to upper case bases
        3. collect the species
        4. collect the number of species, length of maf block
        5. collect species IDs, block_Start coordinates in hg38
        6. IF >= 2 closest relative primate genomes are in the alignment block
            7. parse human-specific gap, human-specific insertions 
            8. convert relative gap/insertion positions to genomic intervals
            9. write the results per block
        10. rezip .maf file. 
    
    """

    ### pre processing ###
    # copy source files to local, get chr number, datapath for output
    CHR, DATA_PATH = copySrcToLocalMaf(INDEX)


    #0
    outfile = os.path.join(DATA_PATH, f"hu-specific-GAPS_INS_{CHR}.bed")

    PRIMATE_GENOMES = ["panTro5", "gorGor5", "rheMac8", "panPan2", "nomLeu3", "ponAbe2"]

    #1
    if os.path.exists(outfile) is False:
        unzipLocalMaf(CHR, DATA_PATH)
        for block in AlignIO.parse(f"{CHR}.maf", "maf"):

            #2 block --> np.array
            align_array = np.array(block)
            align_array = np.char.upper(align_array) # change bases to upper. ignore repeats for now. 

            #3 collect species
            species_vector = []

            #4 count n species, block seq len
            nspecies, maf_block_len = align_array.shape

            #5 collect species, block start coordinates in hg38
            for seq in block:    
                species_vector.append(seq.id.split(".")[0])

                if seq.id == f"hg38.{CHR}":
                    maf_block_start = seq.annotations["start"]

            #6 If there are alignments for 2/4 closest relative primate species 
            n_primates_in_block = len(set(PRIMATE_GENOMES).intersection(set(species_vector)))

            if n_primates_in_block>=2:

                # get primate alignment indexes in sequence block
                primate_block_index = []

                for i, s in enumerate(species_vector):

                    if s in PRIMATE_GENOMES:

                        primate_block_index.append(i)

                #7 parse indels in block (function)
                gaps, ins, sub = parseDelsIns(align_array, species_vector, primate_block_index, maf_block_start)

                lists_to_write = [
                                ("DEL", gaps), 
                                ("INS", ins), 
                                ("SUB", sub)
                                ]

                for mutation_type, mut_list in lists_to_write:

                    if len(mut_list)>0: # if there are bases to write

                        #8 convert gap list to genome coordinate intervals
                        genomic_gaps = getConsecCoorIntervals(mut_list, mutation_type)

                        #9 write the gaps for this block
                        writeInfo(CHR, species_vector, nspecies, maf_block_len, maf_block_start, genomic_gaps, mutation_type, outfile)

        zippery.rezip_file(f"./{CHR}.maf") 
    else:
        print("already wrote", outfile)

if __name__ == "__main__":
    main(sys.argv[1:])