from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import config_readwrite as crw
import os, sys
import numpy as np
import pandas as pd

FASTA = sys.argv[1]
CFN = sys.argv[2]

# functions

def getHandles(fasta):
    """make a dictionary of all the output handles for one fasta file sei run"""
    
    SEI_SRC = "/wynton/home/ahituv/fongsl/bin/sei-framework/"
    SEI_PATH = os.path.join(os.path.split(fasta)[0], "sei_predictions")

    FASTA_CLEAN = os.path.splitext(fasta)[0] + ".clean.fa"
    FASTA_INDEX =  os.path.splitext(FASTA_CLEAN)[0] + ".index.txt"
    
    PADDED  = os.path.splitext(FASTA_CLEAN)[0] + ".sei_padded.fa"
    PATH, HANDLE = os.path.split(PADDED)
    HANDLE = HANDLE.strip(".fa")
    
    
    CHROM_PATH= os.path.join(
        SEI_PATH, "chromatin-profiles-hdf5")
    
    CHROM_OUT = os.path.join(CHROM_PATH, f"{HANDLE}_predictions.h5")
    CLASS_OUT = os.path.join(CHROM_PATH, f"{HANDLE}.raw_sequence_class_scores.npy")
    LABEL_OUT = os.path.join(CHROM_PATH, f"{HANDLE}_row_labels.txt")
    CLASS_TABLE = os.path.join(CHROM_PATH, f"{HANDLE}.raw_sequence_class_scores.table.tsv.gz")
    
    
    path_dict = {
        "FASTA": fasta,
        "FASTA_CLEAN":FASTA_CLEAN, 
        "FASTA_INDEX":FASTA_INDEX,
        "PADDED":PADDED, 
        "PATH": PATH,
        "HANDLE":HANDLE, 
        "SEI_PATH":SEI_PATH, 
        "SEI_SRC":SEI_SRC,
        "CHROM_PATH": CHROM_PATH, 
        "CHROM_PRED": CHROM_OUT,
        "CLASS_PRED" : CLASS_OUT, 
        "LABELS" : LABEL_OUT, 
        "CLASS_TABLE":CLASS_TABLE
    }
    
    return path_dict


def cleanFasta(fasta):
    """
    clean fasta files for sei run by indexing sequence ids
    
    goal 
        index sequence ids 
        remove sequence id duplicates
    
    requires
        Bio.SeqIO.FastaIO import SimpleFastaParser
        getHandles function
        
    input 
        fasta file (str) - path to fasta file
        path_dict (dictionary) - dictionary with files

    method    
        1. get dictionary of file handles
        2. instantiate fasta clean adn index handles to write
        3. open files to write
        4. make empty list to collect unique sequence ids
        5. open fasta and iterate through rows, making uniform sequence ids
            5.1. check to make sure id has not been processed, if not, append to id_list
            5.2. write the index_id and the sequence
            5.3. write the index_id and sequence_id to a separate file.  
        
    Return/writes
        FASTA_CLEAN (str) - fasta file with sequence id indexes
        FASTA_INDEX (str) - .txt file with sequence id index + sequence id
    """
    print("clean and index fasta")
    
    #1
    path_dict = getHandles(fasta)
    
    #2 get files (strs) to write
    FASTA_CLEAN, FASTA_INDEX = path_dict["FASTA_CLEAN"], path_dict["FASTA_INDEX"]
    
    if os.path.exists(FASTA_CLEAN) is False:
        #3 open the write files
        writer, indexer = open(FASTA_CLEAN, "w"), open(FASTA_INDEX, "w")

        #4 list to collect sequence ids
        id_list = []

        #5 open the fasta
        with open(fasta, "r") as handle:
            for i, value in enumerate(SimpleFastaParser(handle)):

                index_id = f"seq.{i}"  # uniform sequence id
                id_, seq = value

                #5.1
                if id_ not in id_list:

                    id_list.append(id_)

                    #5.2.
                    writer.write(f'>{index_id}\n{seq}\n')

                    #5.3.
                    indexer.write(f'{index_id}\t{id_}\n')
                else:
                    print("sequence_id duplicate!", id_)

        writer.close(), indexer.close()
    else:
        print("cleaned fa already", FASTA_CLEAN)
    
    return FASTA_CLEAN, FASTA_INDEX


def seqInsert(seq, insert_start, insert_seq):
    """
    insert sequence fragment (insert_seq) at position (insert_start) within full sequence (seq)
    
    Insert into center if insert_start is None 

    return inserted sequence. 
    """

    insert_size = len(insert_seq)
    
    if insert_start is None:
        
        insert_start = (len(seq)/2) - (insert_size/2)
        #print(insert_start)

    insert_end = insert_start + insert_size  # find center end

    return seq[:int(insert_start)] + insert_seq + seq[int(insert_end):]


def trimSeq(seq, size):
    """ find center of sequence and trim down to size"""
    
    center = len(seq)/2
    start = center - (size/2)
    end = center + (size/2) + 1
    
    return seq[start:end]


def padSeq(fasta):
    
    """ if sequence is shorter than 4096, pad, else trim"""
    
    print("padding/trimming clean fasta to 4096bp")
    
    # get handles
    path_dict = getHandles(fasta)
    
    # set paramters, pad with Ns
    max_len, PAD = 4096, "N"
    
    # out file - padded and clean
    OUT = path_dict["PADDED"]
    
    if os.path.exists(OUT) is False:
        sequences = [s for s in SeqIO.parse(fasta, 'fasta')]

        padded_sequences = []

        for n, seq in enumerate(sequences):
            if len(seq.seq)<max_len:
                padding = PAD*max_len # creating the padding string
                padded_sequences.append(seqInsert(padding, None, seq)) # insert the sequence in the center, append to list
            else:
                padded_sequences.append(trimSeq(seq, max_len))
        SeqIO.write(padded_sequences, OUT, 'fasta')  # write all the sequences
    else:
        print("padded/trimmed already", OUT)

    return OUT


def launchSei(fasta, build, gpu):

    print("launching sei")
    
    path_dict = getHandles(fasta)
    
    SEI_SRC = path_dict["SEI_SRC"]
    SEI_PATH = path_dict["SEI_PATH"]
    
    CHROM_PRED = path_dict["CHROM_PRED"]
    CLASS_PRED = path_dict["CLASS_PRED"]
    
    if os.path.exists(SEI_PATH) is False:
        os.mkdir(SEI_PATH)
        
    print(SEI_PATH, "\n", CLASS_PRED)

    GPU_BOOL = "True" if gpu is True else "False"

    SCRIPT = os.path.join(SEI_SRC, "sarah_scripts/launch_qsub.py")

    cmd = " ".join(['python',
           SCRIPT,
           fasta,
           build,
           SEI_PATH,
           GPU_BOOL,
           CHROM_PRED, 
           ])

    if os.path.exists(CLASS_PRED) is False:
        print(cmd)
        os.system(cmd)
    else:
        print('ran sei already', CLASS_PRED)

def returnSequenceClassLabels():
    file = "/wynton/home/ahituv/fongsl/bin/sei-framework/sequence_class_labels.csv"
    lab = pd.read_csv(file)

    return lab


def processLabel(label_file, index):
    """
    input 
        label_file  (str) - path with the labels for the sequences run through the DNN

    Method
        1. opens sequence label file as pd dataframe
        2. drop the index column
        3. make UCSC genome browswer coordinates
        4. avg pct? 
        5. "bin"?
    """

    # read label file as pd dataframe
    lab = pd.read_csv(label_file, sep='\t')
    lab = lab.drop(columns=["index"])  # redundant index column

    ind = pd.read_csv(index, sep='\t', header=None, names=["name", 'tile.coor'])
    lab = pd.merge(lab, ind)
    return lab[["tile.coor"]]


def makeTable(class_pred, labels, index_file, table_file):

    # get sequence class labels from function above
    seqClass = returnSequenceClassLabels()
    
    # open npy data
    data = np.load(class_pred, allow_pickle=True)

    # turn into pd dataframe
    df = pd.DataFrame(data)

    # rename columns
    df.columns = list(seqClass["#Sequence class label"])[:-1]

    # process labels file
    lab = processLabel(labels, index_file)

    # add labels and data together
    df = pd.merge(lab, df, left_index=True, right_index=True)

    # write table to table file. 
    df.to_csv(table_file, sep='\t', index=False, compression="gzip")
    
    return df
        
def main(argv):
    
    # load config
    config, cfn = crw.read(CFN)
    
    FASTA_CLEAN, FASTA_INDEX = cleanFasta(FASTA)
    
    # sequence padding w n
    PADDED_FASTA = padSeq(FASTA_CLEAN)

    # launch sei
    GPU = True
    launchSei(PADDED_FASTA, "hg38", GPU)  
    
    # write to config
    path_dict = getHandles(FASTA)
    
    SECTION = "sei"
    for key, value in path_dict.items():
        config[SECTION][key] = value

    crw.write(config, cfn)

    """
    # get sequence class labels. See Methods section of Chen 2022 for interpretation of these PC labels.
    # apparently labels >40 are low active/heterochromatin. 
    # Make up <2% of the genome. But 2% of the genome can still be significant.
    """
    if os.path.exists(path_dict["CLASS_PRED"]) is True:
        makeTable(path_dict["CLASS_PRED"], path_dict["LABELS"], FASTA_INDEX, path_dict["CLASS_TABLE"])
    else:
        print("waiting for sei to finish. make table later")
    
if __name__ == "__main__":
    main(sys.argv[1:])    