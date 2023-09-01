import argparse
import os, sys

parser = argparse.ArgumentParser()

parser.add_argument('inputfile', type=str, help=".fa, .vcf, .bed file") 
parser.add_argument('build', type=str, default="hg38", help="hg38, hg19 genome build") 
parser.add_argument('outdir', type=str, default = f"{os.getcwd()}", help="output directory") 

args = parser.parse_args()

FILE, BUILD, OUTDIR,= args.inputfile, args.build, args.outdir

QSUB =  "/wynton/home/ahituv/fongsl/bin/sei-framework/interpret_sei.sh"

cmd = ["qsub", QSUB, FILE, BUILD, OUTDIR]
os.system(" ".join(cmd))