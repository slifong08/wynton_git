
import glob
import os
import subprocess
import sys

# append path
sys.path.append("/wynton/home/ahituv/fongsl/tools/py_/")
import config_readwrite as crw

config_tag = sys.argv[1]

config, n = crw.read_config(config_tag)

# nullomer GENCODE intersections
ANNOT = config["GENCODE"]["ANNOT"]
OVERLAP = config[f"DATAxGENCODE"]["OVERLAP"]  
NOOVERLAP_REF = config[f"DATAxGENCODE"]["NOOVERLAP_REF"]

# nullomer shuffles
EX_EXP_SHUF = config["SHUFFLE"]["ex_exp_star"]
NOEX_EXP_SHUF = config["SHUFFLE"]["noex_exp_star"]


runlist =[
        config["DATAxGENCODE"]["OVERLAP"]
        ]

# ## shuffled exonic, non-exonic; expand; extract bw values

runlist.extend(glob.glob(config["SHUFFLE"]["ex_exp_star"]))
runlist.append(config["DATAxGENCODE"]["NOOVERLAP_REF"])
runlist.extend(glob.glob(config["SHUFFLE"]["noex_exp_star"]))


# write QSUB section. + ARRAY object
crw.check_section(config, "QSUB")

ARRAY = os.path.join(os.getcwd(), f"array-{ANNOT}.txt")

config["QSUB"]["ARRAY"] = ARRAY

crw.write_config(config, n)


# write the array file with numbers. 

with open(ARRAY, "w") as writer:
    for i, file_name in enumerate(runlist):
        j = i+1 # 1-index job? 
        writer.write(f'{j}\t{file_name}\n')
    writer.close()
                     