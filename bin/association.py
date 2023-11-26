#!/usr/bin/env python
# coding: utf-8

# # MPRA flow
# 20230829 
# sarahfong
# 
# 
#     adapted for python + qsub on wynton
#     because nextflow does not play well with wynton
#     
#     Requires environment modules from wynton
#         Loading samtools, bwa, picard
#         
#         ml load CBI picard/2.27.5  samtools/1.18  bwa/0.7.17



import argparse

from Bio.SeqIO.FastaIO import SimpleFastaParser
import glob
import gzip
import os, sys
import numpy as np
import subprocess as sp

sys.path.append("/wynton/group/ahituv/fongsl/tools/py_")
import config_readwrite as crw


# # Args


arg_parser = argparse.ArgumentParser(description="run MPRA flow with python")

arg_parser.add_argument("path", type=str, help='full path to fastq read directory')
arg_parser.add_argument("name", type=str, help='name of experiment')

args = arg_parser.parse_args()

# parse args into variables
PATH = args.path
NAME = args.name


# # DEV PARAMS


# read config
CONFIG_NAME = os.path.join(PATH,"assoc", NAME + ".config.ini")

# read 
config, cfn = crw.read(CONFIG_NAME) 

# new section - params
section = "params"

PATH = config[section]["assoc_path"]
NAME = config[section]["name"]
BASEDIR = config[section]["basedir"] 
DESIGN = config[section]["design"]
LABELS = config[section]["labels"] 


# new section - fastq
section = "fastq"

INS = config[section]["FASTQ_INSERT"]
INSPE = config[section]["FASTQ_INSERTPE"]
FASTQ_BC = config[section]["barcodes"]

print(DESIGN, LABELS)
    
"""
PATH = "/wynton/group/ahituv/fongsl/projects/nullomers/data/20230815_nullomer_MPRA/test"
NAME = "TEST"
DESIGN = os.path.join(os.path.dirname(PATH), "15mer.fo.pam.scaffold.ext200.library.TWIST.fa")
INS = "SF_asso_S1_R1_001.fastq.TEST.gz"
INSPE ="SF_asso_S1_R4_001.fastq.TEST.gz"
FASTQ_BC = "SF_asso_S1_R2_001.fastq.TEST.gz"
LABELS = None


os.chdir(PATH)

### full dataset run

PATH = "/wynton/group/ahituv/fongsl/projects/nullomers/data/20230815_nullomer_MPRA/"
NAME = "nullomer_mpra"
DESIGN = os.path.join(os.path.dirname(PATH), "15mer.fo.pam.scaffold.ext200.library.TWIST.fa")
DESIGN = os.path.join(PATH, "15mer.fo.pam.scaffold.ext200.library.TWIST.fa")
INS = "SF_asso_S1_R1_001.fastq.gz"
INSPE ="SF_asso_S1_R4_001.fastq.gz"
FASTQ_BC = "SF_asso_S1_R2_001.fastq.gz"
LABELS = None


os.chdir(PATH)


"""


# # functions

# # launchCmd


def launchCmd(cmd, **kwargs):

    out = kwargs.get('out', None)
    
    # if cmd is a list and not a string. 
    if type(cmd) is list:
        cmd = " ".join(cmd)

    print("\nRUNNING\n\t", cmd, "\nOUT\n\t", out)
    
    os.system(cmd)
    
    
# ## count_bc


def countBC(fastq_bc, design, labels):
    
    if labels != "None":
        """
        * count fastq and bam length remove the illegal regex characters
        * and make design file
        * contributions: Gracie Gordon & Max Schubach
        */

        """
        cmd = """awk '{gsub(/\\[/,"_")}1' %s > t_new_label.txt""" % labels
        launchCmd(cmd)

        cmd = """awk '{gsub(/\\]/,"_")}1' t_new_label.txt > label_rmIllegalChars.txt"""  #fixed_label
        launchCmd(cmd)

        cmd = """awk '{gsub(/\\[/,"_")}1' %s > t_new_design.txt""" % design
        launchCmd(cmd)

        cmd = """awk '{gsub(/\\]/,"_")}1' t_new_design.txt > design_rmIllegalChars.fa"""  #fixed_design
        launchCmd(cmd)

        cmd = f"zcat {fastq_bc} | wc -l  > count_fastq.txt"  # bc_ch
        launchCmd(cmd)
        
    elif labels == "None":
        """
            if (params.label_file == null) {
            process 'count_bc_nolab' {
                tag 'count'
                label 'shorttime'
                publishDir "${params.outdir}/${params.name}", mode:'copy'

                input:
                    file(fastq_bc) from params.fastq_bc_file
                    file(design) from params.design_file
                output:
                    file 'count_fastq.txt' into bc_ch
                    file "label_rmIllegalChars.txt" into fixed_label
                    file "design_rmIllegalChars.fa" into fixed_design
                shell:

                    #!/bin/bash
                    #CREATE LABEL FILE and remove illegal regex characters
                    awk -F'\t' 'BEGIN {OFS = FS} NR%2==1 {print substr(\$1,2,length(\$1)),"na"}' $design > labels.txt
                    awk '{gsub(/\\[/,"_")}1' labels.txt > t_new_label.txt
                    awk '{gsub(/\\]/,"_")}1' t_new_label.txt > label_rmIllegalChars.txt


                    awk '{gsub(/\\[/,"_")}1' $design | \\
                    awk '{gsub(/\\]/,"_")}1' | \\
                    sed 's/\\r//g' > design_rmIllegalChars.fa

                    zcat $fastq_bc | wc -l  > count_fastq.txt

        """
        new_design_writer = open("design_rmIllegalChars.fa", 'w')
        new_label_writer = open("label_rmIllegalChars.txt", "w")
        
        # use python to open and parse file
        with open(design, "r") as reader:
            for value in SimpleFastaParser(reader):
                label, seq = value
                
                # replace any slashes with underscores... because? 
                # remove open brackets forward slashes close brackets
                new_label = ((label.replace("[", "_")).replace("/", "_")).replace("]", "_")  

                # write
                new_label_writer.write(f"{new_label}\tna\n")
                new_design_writer.write(f">{new_label}\n{seq}")
                
        new_design_writer.close()
        new_label_writer.close()

        cmd = [" sed 's/\\r//g'", "design_rmIllegalChars.fa"]
        if os.path.exists("design_rmIllegalChars.fa") is False:
            launchCmd(cmd, out="design_rmIllegalChars.fa")

        # count fastq lines
        cmd = ["zcat", 
               fastq_bc, 
               "| wc -l  > count_fastq.txt"
              ]
        if os.path.exists("count_fastq.txt") is False:
            launchCmd(cmd, out="count_fastq.txt")



# ## create BWA


def create_BWA_ref():
    """
    /*
    * STEP 1: Align
    * Process 1A: create BWA reference
    * contributions: Gracie Gordon
    */


    process 'create_BWA_ref' {
        tag "make ref"
        label 'shorttime'

        conda 'conf/mpraflow_py36.yml'

        input:
            file(design) from fixed_design
            file(label) from fixed_label
        output:
            file "${design}.fai" into reference_fai
            file "${design}.bwt" into reference_bwt
            file "${design}.sa" into reference_sa
            file "${design}.pac" into reference_pac
            file "${design}.ann" into reference_ann
            file "${design}.amb" into reference_amb
            file "${design}.dict" into reference_dict
        shell:

            #!/bin/bash
            bwa index -a bwtsw $design
            samtools faidx $design
            picard CreateSequenceDictionary REFERENCE=$design OUTPUT=$design".dict"
            """
    section = "countBC"
    DESIGN = config[section]["new_design_noillegal"]
    LABEL = config[section]["new_label_noillegal"]
    cmd = [
        "bwa index -a bwtsw",
        DESIGN,
    ]

    if os.path.exists(DESIGN +".bwt") is False:
        launchCmd(cmd)

    cmd = [
        "samtools faidx",
        DESIGN,
    ]
    if os.path.exists(DESIGN + ".ann") is False:
        launchCmd(cmd)

    cmd = [
        "picard CreateSequenceDictionary REFERENCE=",
        DESIGN,
        "OUTPUT=",
        f"{DESIGN}.dict"
    ]
    if os.path.exists(DESIGN + '.dict') is False:
        launchCmd(cmd, out=(DESIGN + '.dict'))

# ## chunking fastq - can't just split, need to use iterator

def makeSplitFq():
    COUNT_FQ_FN = config["countBC"]["count_fastq"]
    NLINES = config["params"]["split"]
    
    # get the number of fastq lines. 
    with open(COUNT_FQ_FN, "r") as reader:
        for line in reader:
            COUNT_FQ = int(line)
            break
    
    # determine splits from the bc count and desired number of lines
    NSPLITS = int(round(int(COUNT_FQ)/int(NLINES),0))
    print("fastq lines", COUNT_FQ, "desired line split", NLINES, "nsplits", NSPLITS)
    
    fs_ins, fs_inspe = [], []

    for i in np.arange(NSPLITS):
        
        num = int(1+i)
        file_ins, file_inspe= f"-o {num}.group_INS.fq", f"-o {num}.group_INSPE.fq"
        fs_ins.append(file_ins), fs_inspe.append(file_inspe)

    return " ".join(fs_ins), " ".join(fs_inspe)

def chunkingFastQ():
    """/*
    *CHUNKING FASTQ
    */

    Channel
        .fromPath(params.fastq_insert_file)
        .splitFastq( by: params.split, file: true )
        .set{ R1_ch }

    if (params.fastq_insertPE_file != null) {
        Channel
            .fromPath(params.fastq_insertPE_file)
            .splitFastq( by: params.split, file: true )
            .set{ R3_ch }
    }

    /*
    """
    NLINES = config["params"]["split"]
    INS = config["fastq"]["FASTQ_INSERT"]
    INSPE = config["fastq"]["FASTQ_INSERTPE"]
    FASTQ = True     
    
    ### fastq splitter  - if quality control is not needed? 
    ### fastp - when quality control is desired. 
    
    # check that this has not been done. 
    if len(glob.glob("*.group_INS.fq")) ==0:
        if FASTQ is True:
            # output file string
            fs_ins, fs_inspe = makeSplitFq()

            split_pairs = [(INS, fs_ins), 
                           (INSPE, fs_inspe)
                          ]
            
            for insert, split_insert_string in split_pairs:
                cmd_fastq = ["fastqsplitter", 
                                 "-i", 
                                 insert, 
                                 split_insert_string
                                ] 
                launchCmd(cmd_fastq, out="*.group_INS|INSPE.fq")
        else:    
            # path to binary
            FASTP = "/wynton/group/ahituv/bin/fastp"
            
            # command
            cmd_fastp = [
                        FASTP,  
                        "-i",
                        INS,
                        "-I", 
                        INSPE, 
                        "--out1=group_INS.fq", 
                        "--out2=group_INSPE.fq",
                        f"--split_by_lines={NLINES}",
                        "-z 1", 
                        "--disable_quality_filtering", 
                        "--disable_length_filtering", 
                        "--disable_adapter_trimming", 
                        "--dont_eval_duplication", 
                        ]
            launchCmd(cmd_fastp, out="*.group_INS|INSPE.fq")



# ## mergePairedEndReads

 ## mergePairedEndReads

def mergePairedEndReads(fastq_insert, fastq_insertPE):
    """/*
    *Process 1B: merge Paired end reads
    * contributions: Gracie Gordon
    */
    if (params.fastq_insertPE_file != null) {
        process 'PE_merge' {
            tag 'merge'
            label 'shorttime'

            conda 'conf/mpraflow_py36.yml'

            input:
             file(fastq_insert) from R1_ch
            file(fastq_insertPE) from R3_ch
        output:
            file "*merged.fastqjoin" into mergedPE
        shell:
         fastq-join $fastq_insert $fastq_insertPE -o ${fastq_insert}_merged.fastq
    

    had to clone and make fastq-join from ea.utils
    cd /wynton/group/ahituv/fongsl/src/MPRAflow/src
    git clone https://github.com/ExpressionAnalysis/ea-utils.git
    cd  ./ea-utils/clipper
    make install
    """

    print("mergePairedEndReads")
    
    name=fastq_insert.strip("_INS.fastq")

    cmd = [
            #"/wynton/group/ahituv/fongsl/src/MPRAflow/src/ea-utils/clipper/fastq-join",
            "/wynton/group/ahituv/bin/fastq-join",
            fastq_insert,
            fastq_insertPE, 
            "-o" ,
            f"{name}_merged.fastq"
            ]
    if os.path.exists(f"{name}_merged.fastqjoin") is False:
        launchCmd(cmd, out=f"{name}_merged.fastqjoin")

    if os.path.getsize(f"{name}_merged.fastqjoin") > 0:
        cmd = ["gzip",
               fastq_insert,
               fastq_insertPE,
              ]
        
        launchCmd(cmd, out=f"gzip {fastq_insert} {fastq_insertPE}")
        print(" ".join(cmd))
        
    return f"{name}_merged.fastqjoin"


# ## align W BWA  


def alignWithBWA(chunk):
    """
    /*
    * Process 1C: align with BWA
    * contributions: Gracie Gordon
    */

    //paired ends
    if (params.fastq_insertPE_file != null) {
    process 'align_BWA_PE' {
        tag "align"
        label 'longtime'

        conda 'conf/mpraflow_py36.yml'

        input:
            file(design) from fixed_design
            file(chunk) from mergedPE
            val(name) from params.name
            file(reference_fai) from reference_fai
            file reference_bwt from reference_bwt
            file reference_sa from reference_sa
            file reference_pac from reference_pac
            file reference_ann from reference_ann
            file reference_amb from reference_amb
            file reference_dict from reference_dict
        output:
            file "${name}.${chunk}.sorted.bam" into s_bam
            file '*count_bam.txt' into bam_ch
        shell:
    """
    section = "countBC"
    DESIGN = config[section]["new_design_noillegal"]
    
    section = "params"
    NAME = config[section]['name']

    section = "alignWithBWA"
    SBAM = f"{NAME}.{chunk}.sorted.bam"
    CBAM = f'{chunk}.count_bam.txt'
    
    if config["fastq"]["FASTQ_INSERTPE"] != None:
        
        # align chunk to design fasta
        cmd = [
                "bwa mem",
                DESIGN, chunk,
                "| samtools sort -o",
                SBAM,
                ]
        if os.path.exists(SBAM) is False:
            launchCmd(cmd, out=SBAM)

        print('bam made')

        # visual inspection of alignment
        cmd = [
                "samtools view",
                SBAM, "| head"
                ]
        
        #launchCmd(cmd)

        cmd = [
                "samtools view", 
                SBAM, 
                "|",
                f"wc -l >", 
                CBAM
                ]
        if os.path.exists(CBAM) is False:
            
            launchCmd(cmd, out=CBAM)
        
        """
        I DID NOT WRITE THIS IN PYTHON YET
        
        else {
            //single end
            process 'align_BWA_S' {
                tag "align"
                label 'longtime'

            conda 'conf/mpraflow_py36.yml'

            input:
                file(design) from fixed_design
                file(chunk) from R1_ch
                val(name) from params.name
                file(reference_fai) from reference_fai
                file reference_bwt from reference_bwt
                file reference_sa from reference_sa
                file reference_pac from reference_pac
                file reference_ann from reference_ann
                file reference_amb from reference_amb
                file reference_dict from reference_dict
            output:
                file "${name}.${chunk}.sorted.bam" into s_bam
                file '*count_bam.txt' into bam_ch
            shell:

                bwa mem $design $chunk | samtools sort - -o ${name}.${chunk}.sorted.bam
                echo 'bam made'
                samtools view ${name}.${chunk}.sorted.bam | head
                samtools view ${name}.${chunk}.sorted.bam | wc -l > ${chunk}_count_bam.txt


        """

# ##### single insert script is unfinished

# ## collect Fastq chunks


def collectFastqChunks():
    """
    /*
    *COLLCT FASTQ CHUNCKS
    */

    process 'collect_chunks'{
        label 'shorttime'

        conda 'conf/mpraflow_py36.yml'

    input:
        file sbam_listFiles from s_bam.collect()
        file count_bamFiles from bam_ch.collect()
    output:
        file 's_merged.bam' into s_merge
        file 'count_merged.txt' into ch_merge
    script:
        count_bam = count_bamFiles.join(' ')
        sbam_list = sbam_listFiles.join(' ')
    shell:

        #collect sorted bams into one file
        samtools merge all.bam $sbam_list
        samtools sort all.bam -o s_merged.bam

        #collect bam counts into one file

        samtools view s_merged.bam | wc -l > count_merged.txt
    """
    
    section = "alignWithBWA"
    SBAM = config[section]["s_bam"]
    CBAM = config[section]["c_bam"]  # aka bam.ch
    SBAM_MERGE = config[section]["s_bam_merged"]
    CBAM_MERGE = config[section]["c_bam_merged"]
    
    # collect all sbam and cbam files as a string from a list. 
    sbam_listFiles = " ".join(glob.glob(SBAM))
    cbam_listFiles = " ".join(glob.glob(CBAM))
    
    #collect sorted bams into one file
    cmd = [
            "samtools merge all.bam",
            sbam_listFiles
            ]
    if os.path.exists("all.bam") is False:
    
        launchCmd(cmd, out="all.bam")
     
    cmd = [
            "samtools sort all.bam -o",
            SBAM_MERGE
            ]
    if os.path.exists(SBAM_MERGE) is False:
        launchCmd(cmd, out=SBAM_MERGE)
        
    #collect bam counts into one file

    cmd = [
           "samtools view",
           SBAM_MERGE,
           "| wc -l >",
           CBAM_MERGE
          ]
    if os.path.exists(CBAM_MERGE) is False:
        launchCmd(cmd, out=CBAM_MERGE)

# ## assign barcodes


def assignBarcodes():
    """
    /*
    * Assign barcodes to element sequences
    * contributions: Sean Whalen
    */

    process 'map_element_barcodes' {
        tag "assign"
        label "shorttime"
        publishDir "${params.outdir}/${params.name}", mode:'copy'

        conda 'conf/mpraflow_py36.yml'

        input:
            val(name) from params.name
            val(mapq) from params.mapq
            val(baseq) from params.baseq
            val(cigar) from params.cigar
            file(fastq_bc) from params.fastq_bc_file
            file count_fastq from bc_ch
            file count_bam from ch_merge
            file bam from s_merge
        output:
            file "${name}_coords_to_barcodes.pickle" into map_ch
            file "${name}_barcodes_per_candidate-no_repeats-no_jackpots.feather" into count_table_ch
            file "${name}_barcode_counts.pickle"
        shell:

            echo "test assign inputs"
            echo ${mapq}
            echo ${baseq}
            echo $fastq_bc
            zcat $fastq_bc | head

            echo ${count_fastq}
            echo ${count_bam}
            cat ${count_fastq}
            cat ${count_bam}

            python ${"$baseDir"}/src/nf_ori_map_barcodes.py ${"$baseDir"} ${fastq_bc} ${count_fastq} \
            $bam ${count_bam} ${name} ${mapq} ${baseq} ${cigar}
    """

    print("echo test assign inputs")
    
    section = "params"
    NAME =config[section]["name"] 
    BASEQ=config[section]["baseq"]
    MAPQ=config[section]["mapq"]
    CIGAR = config[section]["cigar"] 
    BASEDIR = config[section]["basedir"]
    
    section="fastq"
    FASTQ_BC = config[section]["barcodes"]
    
    section = "countBC"
    COUNT_FASTQ = config[section]["count_fastq"]
    
    section = "alignWithBWA"
    SBAM_MERGE = config[section]["s_bam_merged"]
    CBAM_MERGE = config[section]["c_bam_merged"]
    
    COUNT_BAM = CBAM_MERGE
    BAM = SBAM_MERGE
    
    print("\nbaseq", BASEQ)
    print("\nmapq", MAPQ)

    print("\nfastq_bc", FASTQ_BC)
    
    cmd = [
           "\nzcat", 
           FASTQ_BC,
           "| head"
          ]
    #launchCmd(cmd)
    
    print("Count fastq", COUNT_FASTQ)
    cmd = ["cat", COUNT_FASTQ]
    launchCmd(cmd)
    
    print("count_bam", COUNT_BAM)
    cmd = ["cat", COUNT_BAM]
    launchCmd(cmd)


    cmd = [
            f"python {BASEDIR}/src/nf_ori_map_barcodes.py" , 
            BASEDIR, FASTQ_BC, COUNT_FASTQ, 
            BAM, COUNT_BAM, NAME, MAPQ, BASEQ, CIGAR
            ]
    launchCmd(cmd)

# ## filter barcodes

def filterBarcodes():
    """
    /*
    * Filter barcodes for minimum coverage and unique mapping
    * contributions: Gracie Gordon
    */

    process 'filter_barcodes' {
        tag "$filter"
        label "shorttime"
        publishDir "${params.outdir}/${params.name}", mode:'copy'

        conda 'conf/mpraflow_py36.yml'

        input:
            val(min_cov) from params.min_cov
            val(min_frac) from params.min_frac
            val(out) from params.name
            file(map) from map_ch
            file(table) from count_table_ch
            file(label) from fixed_label
        output:
            file "${out}_filtered_coords_to_barcodes.pickle"
            file "${out}_original_counts.png"
            file "original_count_summary.txt"
            file "${out}_filtered_counts.png"
            file "filtered_count_summary.txt"

        shell:

            python ${"$baseDir"}/src/nf_filter_barcodes.py ${out} ${map} ${table} \
            ${min_cov} ${min_frac} $label
    """
    

    # section - params
    section = "params"
    BASEDIR = config[section]["basedir"]
    NAME = config[section]["name"] 
    MIN_COV = config[section]["min_cov"]
    MIN_FRAC = config[section]["min_frac"]
    
    # relabel
    OUT = NAME
    
    # section - assignBarcodes
    section = "assignBarcodes"
    MAP = config[section]['bc_pickle_coor'] 
    TABLE = config[section]['bc_clean'] 
    config[section]['bc_pickle_count'] = f"{NAME}_barcode_counts.pickle"

    # section - countBC
    section = "countBC"
    LABEL = config[section]["new_label_noillegal"]
    
    cmd =  [f"python {BASEDIR}/src/nf_filter_barcodes.py", 
            OUT, MAP, TABLE, 
            MIN_COV, MIN_FRAC, LABEL
           ]
    launchCmd(cmd)


# # Main


def main(argv):

    os.chdir(PATH)

    countBC(FASTQ_BC, DESIGN, LABELS)

    create_BWA_ref()

    chunkingFastQ()

    for n in (glob.glob("*.group_INS.fq")):
        num = (os.path.split(n)[1]).split(".group")[0]
        chunk_ins_n, chunk_inspe_n = f"{num}.group_INS.fq", f"{num}.group_INSPE.fq"

        print(num, chunk_ins_n, chunk_inspe_n)

        merged_chunk = mergePairedEndReads(chunk_ins_n, chunk_inspe_n)

        alignWithBWA(merged_chunk)


    # glob for chunks? 
    collectFastqChunks()

    assignBarcodes()
    # could take a long time. 

    filterBarcodes()  #also takes a long time

if __name__ == "__main__":
     main(sys.argv[1:])





