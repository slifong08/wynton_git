import os
import subprocess
import sys

# append path
sys.path.append("/wynton/home/ahituv/fongsl/tools/py_/")

# import config reader
import config_readwrite as crw



"""

SET PARAMETERS 

    - DATAPATH (STR) - path to store data
    
    - ANNOTATION (STR) - GENCODE annotation to split nullomers on. 
        - possible GENCODE annotations:

            'exon', 
            'gene', 
            'transcript', 
            'UTR', 
            'start_codon', 
            'CDS', 
            'stop_codon', 
            'Selenocysteine'

    - MUTATION (STR) - input nullomer file (.txt)
    - BUILD (STR) - genome build for these experiments
    
"""


DATA_PATH = "/scratch/fongsl/nullomers/data/"
ANNOTATION = "exon" 
MUTATION = os.path.join(DATA_PATH, "mutations.txt")
BUILD = "hg38"

"""
# step 0 - make config
"""

MAKE_CONFIG = "0_make_config.py"

# pipe function to command line using subprocess

result = subprocess.run(['python', MAKE_CONFIG, ANNOTATION, DATA_PATH, MUTATION, BUILD], stdout=subprocess.PIPE)

# get config name w str split
config =(str(result.stdout).strip("\\n'")).strip("b'")# result.stdout
print(config)

"""
# step 1 - split gencode up by annotation 
"""

SPLIT_GENCODE = "1_extract_gencode.py"

subprocess.run(['python', SPLIT_GENCODE, config], stdout=subprocess.PIPE)
"""
# step 2 - format mutation file into .bed (uniq loci only), .vcf 
"""

FORMAT_MUTATION = "2_format_mutation_file.py"

subprocess.run(['python', FORMAT_MUTATION, config], stdout=subprocess.PIPE)

"""
# step 3 - remove repeatmasker elements
#- 19k mutations after removing repeatmasker
"""

REMOVE_REPEAT_MASKER = "3_remove_repeatmasker.py"

subprocess.run(['python', REMOVE_REPEAT_MASKER, config], stdout=subprocess.PIPE)

"""
# step 4 - separate mutations into coding, non-coding
"""

SEPARATE = "4_separate_annot_mutations.py"

subprocess.run(['python', SEPARATE, config], stdout=subprocess.PIPE)

"""
# step 5 - shuffle

#- create non-coding bkgd exclusion file
"""

SEPARATE = "5_shuffle_annots.py"
ITERS = "500"
subprocess.run(['python', SEPARATE, config, ITERS], stdout=subprocess.PIPE)

"""
# step 6 - phastCons element overlap
"""

PHASTCONS = "6_phastCons.py"

subprocess.run(['python', PHASTCONS, config], stdout=subprocess.PIPE)

"""
# step 7 - Phylop
"""

PHYLOP_EXP = "7-phylop_expand.py"
NBASES_EXP = "500"

subprocess.run(['python', PHYLOP_EXP, config, NBASES_EXP], stdout=subprocess.PIPE)