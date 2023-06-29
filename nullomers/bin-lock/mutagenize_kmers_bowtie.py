"""
Take kmer spectra, nullomer spectra, mutagenize sequences to get kmer spectra of the mutants. 

input
    mer (str) - identity of the kmer length
    nmuts (str) - number of positions to mutagenize simultaneously. 
    cell line (str) - cell ine
    build (str) - genome build
    
"""

import argparse
from Bio.SeqUtils import gc_fraction

from collections import Counter
import glob
import gzip
from itertools import product, combinations
import numpy as np
import os
import sys
import subprocess as sp

sys.path.append("/wynton/home/ahituv/fongsl/tools/py_")
import config_readwrite as crw

parser = argparse.ArgumentParser()
parser.add_argument("cellline", type=str, help = "cellline")
parser.add_argument("kmer_length", type=int, help = "length of kmer")
parser.add_argument("n_mut", type=int, help = "number of mutations to make")
parser.add_argument("build", type=str, help = "genome build")
parser.add_argument("fo_only", type=bool, default=True, help = "keep first order nullomers only?")
parser.add_argument("config", type=str, help = "path to config file")

args = parser.parse_args()
CL, MER, NMUTS, BUILD, FO_ONLY, CONFIG = args.cellline, args.kmer_length, args.n_mut, args.build, args.fo_only, args.config

"""
# dev
CL, MER, NMUTS, BUILD, FO_ONLY = "common", "14", "2", "hs1", True
CONFIG = os.path.join(os.getcwd(), "config.ini")
"""
# functions


def getFirstOrderNulls(kmer_len, build):
    """
    return set of first order nullomers in genome background. 

    input
        kmer_len (int) - length of kmer
        build (str) - genome build of nullomer
    method
        1. make empty set
        2. glob all first order nullomer files
        3. for each file, read and add nullomer sequence to set

    return
        first_order (set) - set of first_order nullomers.
    """
    first_order = set()
    for f in glob.glob(f"/wynton/home/ahituv/fongsl/dna/{build}/kmers/{kmer_len}mers/order*.tsv"):
        with open(f, "r") as reader:
            for line in reader:
                first_order.add(line.strip("\n").split("\t")
                                [0])  # add seq to set

    print(len(first_order))
    return first_order


def getSeqDict(file_list):
    """
    get dictionary of kmer sequences and counts, guides

    input
        file_list (list) - list of files w/ full path

    method 
        1. instantiate collection dict
        2. read file, extract sequence and kmer counts, add to element specific dictionary
    return
        seq_dict (dict) - dictionary of kmers and counts

    """
 
    # print("getting sequence file", file)

    # 1
    seq_dict = {}

    # 2
    for file in file_list:
        key = os.path.split(file)[1]

        if "fwd" in file or "rev" in file:
             key = key.split(".csv.gz")[0]

        else:
            key = (key).split(".")[0]  # get the key AAAA

        if key not in seq_dict:
            seq_dict[key] = {}

        file_dict=seq_dict[key]

        with gzip.open(file, "rt") as reader:

            for line in reader.readlines():
                seq, count = line.strip("\n").split(",")
                file_dict[seq] = count

    return seq_dict


def genKmers(length):
    """
    return all sequence permutations, including repeats (AAAAA, GGGGG, CCCCC etc.)

    require
        itertools.product

    input
        length (int)

    return
        mers_list (list) - list of all nucleotide permutations 
    """
    print("generating kmer space length", length)

    mers = product("ACTG", repeat=length)

    mers_list = list("".join(i) for i in mers)

    return mers_list


def getPosLetterCombos(nmuts, sequence):
    """
    return combinations of (1) indices (2) mutated bases for mutating a sequence

    require 
        itertools

    inputs 
        nmuts (int) - number mutations to make
        sequence (str) - sequence to be mutated

    method
        1. get index combinations based on sequence length and number of mutations to make. 
            This makes a map of all possible combinations of sequences to mutate
            1.1 combinations requires that each index is unique and non-redundant. 
                Order does not matter - 
                    e.g. (2,4,5) is the same as (5,2,4) because indexes 2, 4, and 5, will all be mutated.

        2. get sequence product to mutate at indexes
            2.1 - product allows for repeats of the same base in different positions

    return
        mut_pos (list) - list of positional combinations to mutate
        mut_bases (list) - list of letter combinations to mutate


    """
    # print("making index combinations, nucletide permutations of length", nmuts)

    # 1 index combinations
    mut_pos = list(combinations(np.arange(len(sequence)), nmuts))

    # 2 nucleotide permutations per index combo.
    mut_bases = list(product("ACGT", repeat=nmuts))

    return mut_pos, mut_bases


def buildSeqMut(sequence, mut_pos, mut_bases):
    """
    mutate sequence at position with letter
    multiple positions and letters can be inserted into the sequence simultaneously. 

    input
        sequence (str) - original sequence
        mut_pos (set) - Sets of single, tuple, threeple positional index(es) to mutate
        mut_bases (tuple) - Sets of single, tuple, threeple nucleotide combinations to mutate sequence to

    method
        1. instantiate seqs set to collect mutated sequences, add identity to seq set
        2. per positions to mutate in set
            2.1 per base combinations to mutate at these positions
        3. zip and iterate through positions and bases, mutating input sequence
        4. IF mut_seq != input sequence, then return. Else, skip


    return 
        seqs (set) - set of sequences with mutations 

    """

    # 1
    seqs, mut_seq = set(), ""

    # 2
    for pos in mut_pos:

        # 2.1
        for letters in mut_bases:

            # 3
            for p, l in zip(pos, letters):

                if mut_seq == "":
                    mut_seq = sequence[:p] + l + sequence[p + 1:]

                else:
                    mut_seq = mut_seq[:p] + l + mut_seq[p + 1:]

            # 4
            if mut_seq != sequence:
                seqs.add(mut_seq)
                mut_seq = ""

            else:
                # print('no mut', sequence, mut_seq
                mut_seq = ""
                pass
    # add the identity sequence back into the mix.
    seqs.add(sequence)
    return seqs


def generateMismatchSpectra(sequence, nmuts, kmer_spectra):
    """
    input 
        sequence (str) - sequence to mutate and match with kmer-spectra
        nmuts (int) - max number of mismatches to mutate each kmer sequence by
        kmer_spectra (dict) - dictionary of kmer keys and their frequency count (value)

    require
        getPosLetterCombos function
        buildSeqMut function
        prettifySeq function

    method
        1. instantiate kmer_match, null_match list.

        2. get all combinations of indexes and mutated bases to try. 
            mutate bases will be a single, tuple, or threeple depending on the number of mutations desired.

        3. mutate all positions of the nullomer with that base. 

            3.1 - get mutated sequences using buildSeqMut function
                NOTE - mutBase function will remove identity sequence. See function method step 4.3.

        4. per mutated sequence 
            4.1 look up kmer count of the mutated sequence. Test whether it is a prime. 
            4.2 add any resurfaced nullomer sequences to the resurface dict 
                - that is, resurfaced nullomers are sequences where mutations create nullomers
                - return this dictionary
            4.3 add kmer count to the dictionary as pretty key
        7. Add seq kmer dictionary back into the collection dictionary. 


    return
        kmer_match (list) - all mis-match sequences w/ n muts that are kmers 
        null_match (list) - all mis-match sequences w/ n muts that are NOT kmers (therefor nullomers) 
            every mismatch combination for each sequence 
            and per mismatch combo, the kmer counts for each mismatch

    """ 
    ## print("making kmer mismatch spectra w/ N mutations =", nmuts)
    # 1
    kmer_match, null_match = Counter(), Counter()

    # 2
    mut_pos, mut_bases = getPosLetterCombos(int(nmuts), sequence)  # (1,2) (A,C)

    # 3
    seqs = buildSeqMut(sequence, mut_pos, mut_bases)

    # 4
    for mut_seq in seqs:

        # 4.1 - kmer frequency
        if mut_seq in kmer_spectra[mut_seq[:4]].keys():
            kmer_match[mut_seq] += 1

        # 4.2 - nullomer after mutations, not in kmer-verse
        else:
            null_match[mut_seq] += 1  # append to list


    #print(len(kmer_match), len(null_match), len(seqs))

    return dict(kmer_match), dict(null_match)


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


def buildOtherSeqMuts(original_seq, null_seq):
    """
    mutate original sequence at positions 1 and 2 with two alternative bases, 
    which are not the original bases nor the target nullomer base. These
    are controls for downstream analysis of nullomer mutations
    
    input
        original_sequence (str) - original sequence
        null_sequence (str) - the nullomer sequence
        pos1 (int) - first position in sequence to mutate to other bases
        nucs1 (tuple) - two other nucleotides (as tuple) to mutate position 1 to
        pos2 (int) - second position in sequence to mutate to other bases
        nucs2 (tuple) - two other nucleotides (as tuple) to mutate position 2 to 
        
    method
        1. get all other possible nucleotide combinations that 
        make otherseqs, empty set to collect mutated sequence results. 
        2. iterate through possibilities (4) 
            2.1 set variables for each nucleotide identity at position one and two to be mutated to. 
            n1 - nucleotide identity at position one
                e.g. pos1 = 5, nucs1 = (c,t)
                When n1=0, insert "c" at position 1 (index 5). When n1=1, insert "t" at position 5
            n2 - nucleotide identity at position two
                e.g. pos2 = 9, nucs1 = (a,g)
                When n1=0, insert "a" at position 2 (index 9). When n1=1, insert "g" at position 9
        3. add the mutated sequence to the sequence set
        
    return
        otherseqs (set) - set of other mutated sequences that are not original or nullomer mutations. 
        
    """
    # 1
    combos = []  # combos of other nucleotides. 
    pos = 0 # position counter
    NUC = {"A", "C", "G", "T"}  # set of nucleotides. 

    for o, m in zip(original_seq, null_seq):

        if o != m:
            combos.append((pos, tuple(NUC.difference(set([o,m])))))
        pos +=1  # keep track of the position you're in
    print(combos)

    pos1, nucs1 = combos[0]
    pos2, nucs2 = combos[1]
    
    #1
    otherseqs = set()
    
    #2
    for i in np.arange(4):
        #2.1
        n1 = 0 if i <2 else 1  # toggle first position nucleotide
        n2 = 0 if i%2 ==0 else 1  # toggel second position nucleotide

        #3 build sequence with two mutations to other possible nucleotide combinations
        otherseqs.add(original[:pos1] + nucs1[n1] + original[pos1 + 1:pos2] + nucs2[n2] + original[pos2 + 1:])

    return otherseqs


def makeNullrow(seq, seq_count, null_list, nmuts, coor, datatype):
    """
    return array with information about the original kmer sequence, its frequency, and matching nullomers. 

    inputs
        seq (str) - kmer sequence
        seq_count (int) - count of kmer occurrence in dataset
        null_list (list) -list  of nullomers that are n mutations away
        nmuts (str) - number of mutations between kmer and nullomer
        coor (str) - coordinates of sequence
        datatype (str) -  nullomer, kmer, or kmer_control
        
        

    method
        1. prettify the nullomer mutations (makes it easy to read)
        2. get gc fraction for kmer, nullomers
        3. construct a list of row items 
        
        rownames = ["kmer", 'nKmer', "gc_frac", "relatedNulls", 
                    "nNulls", gc_fracNulls", prettyNulls", "nMuts"
                    ]
        

    return 
        (str) - tab-joined row list + new line
    """

    # 1
    rows = []
    
    gc_seq = round(gc_fraction(seq), 2)

    for n, mut in enumerate(null_list):
        
        pretty_seq = str(prettifySeq(seq, mut))

        gc_mut = str(round(gc_fraction(mut), 2)) 

        rowlist = [seq,  # original kmer
               str(seq_count),  # number of times kmer occurs in set
               str(gc_seq),  # kmer gc fraction
               mut,  # related nullomer sequences
               str(len(null_list)),  # number of related nullomers
               gc_mut,  # gc fraction of related nullomers
               pretty_seq,  # gc fraction of related nullomers
               f"{datatype}-{n}", # to map to fa
                   coor
               ]
        
        rows.append("\t".join(rowlist) + "\n")

    return rows


# WRITE FA files with this process 

def writeFa(seq, null_list, seq_id, datatype):
    
    rows = []
    for n, null_seq in enumerate(null_list):
        rows.append(f">{seq}.{seq_id}.{datatype}-{n}\n{null_seq}\n")
    
    return "".join(rows)

def makeOutfiles(element, outdir, mer, nmuts, fo_only):
    
    if fo_only is True:
        tsv = os.path.join(outdir, f"{element}.{mer}mers.{nmuts}mut.nulls.fo.tsv")
        fa = os.path.join(outdir, f"{element}.{mer}mers.{nmuts}mut.nulls.fo.fa")
        
    else:
        tsv = os.path.join(outdir, f"{element}.{mer}mers.{nmuts}mut.nulls.tsv")
        fa = os.path.join(outdir, f"{element}.{mer}mers.{nmuts}mut.nulls.fa")

    return tsv, fa

def writer(row, outfile):
    with open(outfile, "a") as writer:
        writer.write(row)
        writer.close()
        
def getCoor(key, mer):

    #key, mer= "AACTA", '14'
    COOR_REF = f"/wynton/home/ahituv/fongsl/MPRA/agarwal_2023/kmers/{mer}mers/coor/{key}_coor.csv.gz"
    coor_dict = {}
    with gzip.open(COOR_REF, 'rt') as reader:
        for line in reader:
            parts = line.strip('\n').split(',')
            value_seq, coors = parts[0], parts[1:]

            stripped = set([i.split("_")[0] for i in coors if i != ""]) # take set of coordinates, remove empty, remove +/-
            coor_dict[key+value_seq] = ",".join(stripped)
    return coor_dict

# main


def main(argv):

    """
    mutagenize kmers by n mutations and evaluate if they are still kmers or if they become nullomers

    1. get kmer search space
        1.1 change first order only flag to False if length of kmer is < 14. No first order nullomers below 14bp. 
        
    2. get genome search space, outdir to write to 
        2.1 get first order nullomer set
        
    3. get kmer files, make them into dictionary. 
        add_key = True means file is separated by key sequence (in filename) 
        and value sequences (all sequences sharing that key). Full kmer sequences = key + value sequences. 

    4. get genome files, make those kmers into dictionary. 
        4.1 filter all the nullomer sequences out of the file list
        4.2 add_key = False. Full kmer sequence (with key and value) are already written into genome file. 
        
    5. retrieve outfiles to write (.tsv, .fa)
        5.1 - write to config file. 
        
    6. mutate kmers to identify related nullomers

    7. write nullomers
            7.1. get nullomers that match kmer
            7.2. intersect first order w/ nullomers, check if any are first order. 
            7.3 if nulls (or first order), write .tsv, .fa

    """

    #1
    if CL =="common":
        KMER_REF = f"/wynton/home/ahituv/fongsl/MPRA/agarwal_2023/kmers/{MER}mers/keys/*.csv.gz"
        COOR_REF = f"/wynton/home/ahituv/fongsl/MPRA/agarwal_2023/kmers/{MER}mers/coor/*.csv.gz"
    else:
        KMER_REF = f"/wynton/home/ahituv/fongsl/nullomers/data/lock/{CL}/kmers/{MER}mers/*.csv.gz"

    #1.1
    if int(MER)<14:
        FO_ONLY=False
    else:
        FO_ONLY = True

    #2
    GENOME_REF = f"/wynton/home/ahituv/fongsl/dna/{BUILD}/kmers/{MER}mers/*.csv.gz"
    OUTDIR = f"/wynton/home/ahituv/fongsl/nullomers/data/lock/{CL}/"
    
    #2.1
    FIRST_ORDER = getFirstOrderNulls(MER, BUILD)
    
    #3 open up and get kmer sequence dictionaries.
    KMERS = glob.glob(KMER_REF)

    kmer_seqs = getSeqDict(KMERS)

    print("kmer seqs", list(kmer_seqs.keys())[:3])

    #4 get all kmer sequences in the genome.
    genome_files = glob.glob(GENOME_REF)
    GKMERS = []
    
    #4.1 filter out nullomer.csv files
    [GKMERS.append(k) if "nullomer" not in k else next for k in genome_files]
    
    #4.2
    genome_seqs = getSeqDict(GKMERS)

    print("genome_seqs", list(genome_seqs.keys())[:3])
    

    # 5.1 write to config
    config, cfn = crw.read(CONFIG)
    section = f"{CL}.{MER}mer.{NMUTS}mut"
    crw.check(config, section)

    config[section]["KMER_REF"] = KMER_REF
    config[section]["GENOME_REF"] = GENOME_REF
    config[section]["PATH"] = OUTDIR
   
    crw.write(config, cfn)
    
    # 5 make the out file
    tsv, fa = makeOutfiles(CL, OUTDIR, MER, NMUTS, FO_ONLY)  # make chr-based file

    # 6 mutate kmers to find nullomer sequences.
    for element, seq_dict in kmer_seqs.items():  # element = 'chr1.24114961-24115161.fwd'

        seq_id = element.split(".csv.gz")[0]

        for seq, seq_count in seq_dict.items():  #seq = "ACCGGTTCCACGGT"

            sequence = seq_id + seq # seq_id = first letters of sequence, seq = last letters of sequence

            # mutate one sequence
            related_kmers, related_nullomers = generateMismatchSpectra(
                sequence, int(NMUTS), genome_seqs)

            # 7 if nullomers appear after mutagenesis

            # 7.1
            nulls = [s for s in related_nullomers.keys()]
            control_kmers = [k for k in related_kmers.keys()]

            # 7.2
            first_order = FIRST_ORDER.intersection(set(nulls))

            # First order filter - write ALL nullomers or ONLY first order (fo)?
            #if FO_ONLY is True:
            null_list = list(first_order)
            #else:
                #null_list = nulls

            # 7.3
            if len(null_list) > 0:
                print("write", sequence)
                
                infos = [
                    ([sequence], "kmer"),
                    (null_list, "null"),
                    (list(np.random.choice(control_kmers, size = 30)), 'ctrl')
                    ]

                # get coordinate info
                coor = getCoor(seq_id, MER)
                seq_coor = coor[sequence]

                for seq_list, datatype in infos:
                    farow = writeFa(sequence, seq_list, seq_coor, datatype)
                    rows = makeNullrow(sequence, seq_count, seq_list, NMUTS, seq_coor, datatype)

                    # open writer files
                    fwriter, twriter = open(fa, "a"), open(tsv, "a")

                    # write fa
                    fwriter.write(farow)

                    # write tsv
                    for row in rows:
                        twriter.write(row)

                    twriter.close(), fwriter.close()


if __name__ == "__main__":
    main(sys.argv[1:])