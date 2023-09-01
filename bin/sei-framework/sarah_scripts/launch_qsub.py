import argparse
import os, sys

parser = argparse.ArgumentParser()

parser.add_argument('inputfile', type=str, help=".fa, .vcf, .bed file") 
parser.add_argument('build', type=str, default="hg38", help="hg38, hg19 genome build") 
parser.add_argument('outdir', type=str, default = f"{os.getcwd()}", help="output directory") 
parser.add_argument('gpu', type=str, default = False, help="run GPU?") 
parser.add_argument('inputfile2', type=str, help="chromatin .hdf5 predictions generated from step 1 of sei")


args = parser.parse_args()

FILE, BUILD, OUTDIR, GPU = args.inputfile, args.build, args.outdir, args.gpu,
FILE2 =  args.inputfile2

OUTDIR2 = os.path.join(OUTDIR, "chromatin-profiles-hdf5")

if GPU =="False":
    QSUB =  "/wynton/home/ahituv/fongsl/bin/sei-framework/sarah_scripts/interpret_sei.sh"
    cmd = ["qsub", QSUB, FILE, BUILD, OUTDIR, FILE2, OUTDIR2]
else:
    QSUB =  "/wynton/home/ahituv/fongsl/bin/sei-framework/sarah_scripts/interpret_sei-gpu.sh"
    cmd = ["qsub -q gpu.q", QSUB, FILE, BUILD, OUTDIR, FILE2, OUTDIR2]


os.system(" ".join(cmd))