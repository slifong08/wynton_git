import argparse

import glob
import os, sys
import pandas as pd

sys.path.append("/wynton/home/ahituv/fongsl/tools/py_")
import config_readwrite as crw

def runFimo(fa, jaspar, results_dir):
    # run fimo

    outdir=os.path.dirname(results_dir)
    print(outdir)

    os.chdir(outdir)

    cmd = [
        "fimo",
        jaspar,
        fa, "&& mv" ,
        os.path.join(outdir, "fimo_out"),
        results_dir

    ]
    print(" ".join(cmd))
    os.system(" ".join(cmd))
    
def pamToFa(pam_file):
    # write PAM FA

    PAM_FA =pam_file.strip(".txt") + ".fa"

    if os.path.exists(PAM_FA) is False:
        writer = open(PAM_FA, "w")
        with open(PAM, "r") as reader:
            for i, line in enumerate(reader):
                line=line.strip("\n")
                row=f">{i}\n{line}\n"

                writer.write(row)
        writer.close(), reader.close()
    return PAM_FA
    
    
# write to config
CL = "14mer.firstorder.pam"
config, cfn = crw.read(os.path.join(os.getcwd(), "config.ini"))
section = f"{CL}.FIMO"
crw.check(config, section)
crw.check(config, CL)
PAM="/wynton/home/ahituv/fongsl/dna/hs1/kmers/14mers/14mer.firstorder.pam.purine.nohomopoly.GC.txt"
PAM_FA=pamToFa(PAM)

SEEDPAM="/wynton/home/ahituv/fongsl/dna/hs1/kmers/14mers/SEED_14mer.firstorder.pam.purine.nohomopoly.GC.txt"
SEED_PAM_FA=pamToFa(SEEDPAM)

JASPAR="/wynton/group/ahituv/tfbs_motif/jaspar/JASPAR2022_CORE_non-redundant_pfms_meme.txt"
OUTDIR = "/wynton/home/ahituv/fongsl/nullomers/data/lock/14.fo.pam"
config[section]["jaspar"] = JASPAR
config[CL]["fa"] = PAM_FA
config[CL]["txt"] = PAM
config[CL]["path"]=OUTDIR

FIMO_RESULT_DIR = os.path.join(OUTDIR, f"fimo-{CL}")
FIMO_RESULT = os.path.join(OUTDIR, f"fimo-{CL}", "fimo.tsv")
config[section][f"fimo_dir"] = FIMO_RESULT_DIR
config[section][f"fimo_results"] = FIMO_RESULT

crw.write(config, cfn)

runFimo(SEED_PAM_FA, JASPAR, FIMO_RESULT_DIR)