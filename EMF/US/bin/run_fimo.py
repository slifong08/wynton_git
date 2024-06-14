import os, sys
import pandas as pd

# import args
OUTDIR=sys.argv[1] #./fimo/ALL.ALL
FA=sys.argv[2]
MEME_INDEX=sys.argv[3]

MEME_INDEX_ZERO = int(MEME_INDEX)-1 # minus 1 to account for zero indexing in pandas, one indexing w sge

SRC="/wynton/group/ahituv/bin/meme-5.5.5/src/fimo"

# args for fimo run
FLAGS = " ".join([
                '--best-site',
                '--oc'
                 ])


# get meme file
MEME_PATH="/wynton/group/ahituv/tfbs_motif/jaspar/"
MEME_ARRAY = os.path.join(MEME_PATH, "meme_array.csv")

# lookup meme file from array. 
meme = pd.read_csv(MEME_ARRAY, names=["meme_file"])
meme_fn = meme.iloc[MEME_INDEX_ZERO]["meme_file"]
meme_file = os.path.join(MEME_PATH, meme.iloc[MEME_INDEX_ZERO]["meme_file"])

print(MEME_INDEX_ZERO, meme_file)

outfile = os.path.join(OUTDIR, "US." + meme_fn.strip(".meme") + ".tsv")
print("\n\n", outfile, "\n\n")


# build fimo command
cmd = " ".join([
                SRC, 
                FLAGS,
                OUTDIR, 
                meme_file, 
                FA, 
                ">",
                outfile
                ])

# launch
os.system(cmd)