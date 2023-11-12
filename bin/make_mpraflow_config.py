import argparse
import os, sys
sys.path.append("/wynton/group/ahituv/fongsl/tools/py_")
import config_readwrite as crw

arg_parser = argparse.ArgumentParser(description="run MPRA flow with python")

arg_parser.add_argument("path", type=str, help='full path to fastq read directory')
arg_parser.add_argument("name", type=str, help='name of experiment')
arg_parser.add_argument("INS", type=str, help='Fastq insert reads')
arg_parser.add_argument("PE",  type=str, help='Fastq Paired end reads')
arg_parser.add_argument("barcode", type=str, help='barcodes as demultiplexed i7 fastq reads')
arg_parser.add_argument("design", type=str, help='fasta file of library design')
arg_parser.add_argument("--labels", type=str, required=False,
                        default="None",  help='labels for sequences')


args = arg_parser.parse_args()

# parse args into variables
PATH = args.path
NAME = args.name
DESIGN = args.design
INS, INSPE = args.INS, args.PE
LABELS = args.labels
FASTQ_BC = args.barcode

LABELS = None  # stop gap
ASSOC_PATH = os.path.join(PATH, "assoc")

CONFIG_NAME = os.path.join(ASSOC_PATH, NAME + ".config.ini")

print("\nCONFIG_NAME:", CONFIG_NAME)

# touch
if os.path.exists(CONFIG_NAME) is False:
    
    # make dir
    if os.path.exists(ASSOC_PATH) is False:
        os.mkdir(ASSOC_PATH)
    cmd = ["touch", CONFIG_NAME]
    os.system(" ".join(cmd))
else:
    print("made config already", CONFIG_NAME)

# specify basedir
BASEDIR = "/wynton/group/ahituv/MPRAflow"

# read 
config, cfn = crw.read(CONFIG_NAME) 

# new section - params
section = "params"
crw.check(config, section)

config[section]["path"] = PATH
config[section]["assoc_path"] = ASSOC_PATH
config[section]["name"] = NAME
config[section]["basedir"] = BASEDIR
config[section]["design"] = DESIGN
config[section]["labels"] = str(LABELS)
config[section]["min_cov"]="3"
config[section]["min_frac"]="0.7"
config[section]["baseq"]="30"
config[section]["mapq"]="5"
config[section]["cigar"]="n"
config[section]["split"]="2000000"

# new section - fastq
section = "fastq"
crw.check(config, section)

config[section]["FASTQ_INSERT"] = INS
config[section]["FASTQ_INSERTPE"] = INSPE
config[section]["barcodes"] = FASTQ_BC

# new section - countBC

section = "countBC"
crw.check(config, section)

config[section]["path"]=ASSOC_PATH
config[section]["new_label"] = f"%(path)s/t_new_label.txt"
config[section]["new_label_noillegal"] = f"%(path)s/label_rmIllegalChars.txt"
config[section]["new_design_noillegal"] = f"%(path)s/design_rmIllegalChars.fa"
config[section]["count_fastq"] = f"%(path)s/count_fastq.txt"

# new setion - create BWA.ref
section="createBWA.Ref"
crw.check(config, section)

config[section]["path"]=ASSOC_PATH
config[section]["reference_fai"]= f"%(path)s/{DESIGN}.fai"
config[section]["reference_bwt"]=f"%(path)s/{DESIGN}.bwt"
config[section]["reference_sa"]=f"%(path)s/{DESIGN}.sa"
config[section]["reference_pac"]=f"%(path)s/{DESIGN}.pac"
config[section]["reference_ann"]=f"%(path)s/{DESIGN}.ann"
config[section]["reference_amb"]=f"%(path)s/{DESIGN}.amb"
config[section]["reference_dict"]=f"%(path)s/{DESIGN}.dict"

section = "mergePairedEndReads"
crw.check(config, section)
config[section]["path"]=ASSOC_PATH
config[section]["chunk1"] = f"%(path)s/group_*_INS.fastq"
config[section]["chunk2"] = f"%(path)s/group_*_INSPE.fastq"
config[section]["merged"] = f"%(path)s/group_*_merged.fastqjoin"

section = "alignWithBWA"
crw.check(config, section)
config[section]["path"]=ASSOC_PATH
config[section]["s_bam"]= f"%(path)s/{NAME}.*.sorted.bam"
config[section]["c_bam"] = f"%(path)s/*count_bam.txt"
config[section]["s_bam_merged"]= f"%(path)s/s_merged.bam"
config[section]["c_bam_merged"] = f"%(path)s/count_merged.txt"

section = "assignBarcodes"
crw.check(config, section)

config[section]["path"]=ASSOC_PATH
config[section]['bc_pickle_coor'] = f"%(path)s/{NAME}_coords_to_barcodes.pickle"
config[section]['bc_clean'] = f"%(path)s/{NAME}_barcodes_per_candidate-no_repeats-no_jackpots.feather"
config[section]['bc_pickle_count'] = f"%(path)s/{NAME}_barcode_counts.pickle"
config[section]["all"]=f"%(path)s/{NAME}_barcodes_per_candidate.feather" 
config[section]["no_jackpot"]=f"%(path)s/{NAME}_barcodes_per_candidate-no_jackpots.feather"
config[section]["no_repeats"]=f"%(path)s/{NAME}_barcodes_per_candidate-no_repeats.feather"


section = "filterBarcodes"
crw.check(config, section)

config[section]["path"]=ASSOC_PATH
config[section]["org.counts"]= f"%(path)s/{NAME}_original_counts.png"
config[section]["org.counts.sum"]= f"%(path)s/original_count_summary.txt"
config[section]["fltrd.counts"]=f"%(path)s/{NAME}_filtered_counts.png"
config[section]["fltrd.counts.sum"]=f"%(path)s/filtered_count_summary.txt"

crw.write(config, cfn)