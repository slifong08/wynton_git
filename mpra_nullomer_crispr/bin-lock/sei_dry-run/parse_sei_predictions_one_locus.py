import argparse
import glob
import h5py
import os, sys
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

sys.path.append("/wynton/home/ahituv/fongsl/tools/py_")

import config_readwrite as crw

import plot_params as pp
pp.fonts()


# parse args
parser = argparse.ArgumentParser()
parser.add_argument("file", type=str, help="full path to hdf5 file")
parser.add_argument("cellline", type=str, help = "cellline")
parser.add_argument("kmer_length", type=int, help = "length of kmer")
parser.add_argument("n_mut", type=int, help = "number of mutations to make")

args = parser.parse_args()
FILE, CL, KMER_LEN, NMUTS = args.file, args.cellline, args.kmer_length, args.n_mut


"""
CL, MER, NMUTS = "common", 15, 2
COOR="chr10:133080140-133080340"
PATH = "/wynton/home/ahituv/fongsl/nullomers/data/lock/common/concatmers/sei_predictions/chromatin-profiles-hdf5"
FILE = os.path.join(
    PATH, f"common.{MER}mers.{NMUTS}mut.nulls.fo.pam.CONCAT-10.ext4096.{COOR}_predictions.h5")
"""
###
# FUNCTIONS
###

def readHdF5(filename, path):
    """
    read hdf5 file, dataset
    
    input
        filename (str) - name of hdf5 results
        path (str) - abs path to directory
        
    require 
        h5py
        
    method 
        1. assemble file
        2. read file
        3. get dataset

    return
        file, data
    
    """
    #1
    file = os.path.join(path, filename) 
    
    #2
    f = h5py.File(file, 'r')
    
    #3
    dset = f['data']

    return dset

def getIndexNames():
    
    """
    return list of 21097 index names corresponding to sei prediction features
    
    input
        none
    method 
        1. make a list for the index names
        2. open the index names file 
        3. append name to index name list
    
    return
        idxnames (list) - list of index names
    """
    
    #1
    idxnames = []
    
    #2
    IDXNAMES = "/wynton/home/ahituv/fongsl/bin/sei-framework/seq_prediction_columns.txt"
    with open(IDXNAMES, "r") as reader:
        for line in reader:
            #3
            idxnames.append(line.strip("\n"))
            
    return idxnames

def arrayToDF(dset):
    
    """
    convert np array to pd dataframe 

    input
        dset - np array of arrays
    method
        1. np.vstack dset, turn into dataframe
        2. add column names as index names
   
   return
       df (pd dataframe) - dataframe of the np arrays w reset index
    """
    
    df = pd.DataFrame(np.vstack(dset))
    print(df.shape)
    df.columns = getIndexNames() # function get column names
    
    return df.reset_index()

def getRowNames(row_textfile, path):
    """
    return list corresponding to sample names

    input
        row_textfile (str) - file to read w/ rownames
        path (str) - path to dir w/ text file
    
    method
        1. make empty dictionary to collect rows, instantiate the file
        2. parse through file, add row index, name to dictionary 
        3. make dictionary into pd dataframe
        4. change datatype of index col to int
    
    return 
        df (pd dataframe) - dataframe of row labels + 
            corresponding indexes in prediction file
    """
    
     #1
    rownames = {}
    file = os.path.join(path, row_textfile) 
    
    #2
    with open(file, "r") as reader:
        for line in reader:
             
            idx, name = line.strip("\n").split("\t")
            if idx !="index":
                rownames[idx]=name
    #3            
    df = pd.DataFrame(rownames.items(), columns= ["index_col", "id"])

    #4
    df["index_col"]=df["index_col"].astype(int)
    
    return df


def formatDf(df, names, coor):
    """
    format sei prediction dataframe, 
    
    add 
        row names
        coordinates
        pair_id 
        nullomer 
    fix
        first "id" row, which corresponds to the endogenous scaffold sequence
        
    indput
        df (pd dataframe) - dataframe of sei predictions
        names (pd dataframe) - dataframe of row names and indexes
        coor (str) - coordinates of scaffold sequence
    
    method
        1. merge row names and dataframe predictiosn
        2. drop the index columns. we don't need these anymore. 
        3. fix the id of the first row - endogenous sequence
        4. add columns
            coordinates
            pair_id
            nullomer
    return 
        df (pd dataframe) - formatted dataframe. 
    
    """

    #1 add to dataframe
    df = pd.merge(names, df, how="left", left_on ="index_col", right_on = "index")

    #2 drop the index_col name (redundant)
    df = df.drop(columns=["index_col", "index"])

    # format dataframe
    #3 first item needs to be renamed as the endogenous sequence. 
    df.loc[df.index == 0, "id"] = f"{coor}_-1_concat-endog"
    
    # 4
    df["coor"] = df["id"].apply(lambda x: x.split("_")[0])
    df["pair_id"] = df["id"].apply(lambda x: "_".join((x.split("_")[1:])))

    # add nullomer bool
    df["nullomer"] = False 
    df.loc[df["id"].str.contains("null"), "nullomer" ] = True
    
    return df


def filterCL(df):
    """
    filter dataframe for cell-type specific features
    
    input
        df (pd dataframe) - full dataframe of 21907 predictions + meta columns, rows corresponding to each sequence prediction
        
    method
        1. determin meta-data columns to keep in list. I will add the prediction columns to this list. 
        2. determine cell lines to filter for
        3. per cell line, iterate through 21907 dataframe column names, 
            adding any column name that contains cell line. 
            Uses str search. 
        4. subset dataframe for cl_cols. 
        
    return
        df (pd dataframe) - dataframe reduced to only cellline-specific feature predictions. 
    
    """
    #1
    cl_cols = ["id", "coor","pair_id", "nullomer"]
    
    #2
    CLS = ["K562", "HepG2", " WTC11"]
    
    #3
    for cl in CLS:
        for col in list(df):
            if cl in col:
                cl_cols.append(col)
    
    #4
    df = df[cl_cols].drop_duplicates()
    
    return df

def computeFeatureMean(feature, df):
    """
    compute the mean of the feature for kmers, nullomers, and endogenous sequence
    
    input
        feature (str) - name of feature (e.g. H3K4me3) 
        df (pd dataframe) - prediction datafrme
        
    method
        1. pivot longform dataframe, keeping only feature. 
            Index is each sample, 
            Columns are the tracks measuring that feature 
            (there can be more than one track that measures a feature, like 60-something DHS tracks for k562)
            
            Values will be the predicted probabilities

        2. if there are predictions
        3. per column, calculate the mean of each sequence across multiple track measures for one feature.  
            (except the first, which is the sample id col) 
        4. add feature column
        5. compute difference between each sample pred and endogenous sequence prediction means
        6. compute the rank of each sequence mean 
        7. compute standard score for distribution of sequence means
        
    return 
        dif (pd dataframe) - dataframe of summarized mean scores for the feature, 
            as well as the difference from the endogenous sequence
            
    """
    # pivot dataframe related to feature
    dif = pd.pivot(df.loc[df["feature"]==feature], index="pair_id", columns="track",
                                 values="pred_prob").reset_index()
    if len(dif) >0:

        # calculate the mean score
        dif["mean"] = dif[dif.columns[1:]].mean(axis=1)

        # add feature as a column
        dif["feature"] = feature

        # calculate distance from kmer to each mutation
        dif["endog_dif"] = dif["mean"] - dif.loc[dif["pair_id"].str.contains('endog'), "mean"].iloc[0]

        # rank the means. 
        dif["mean_rank"]=dif["mean"].rank()

        #dif["mean_z"]=stats.zscore(dif["mean"])

        return dif[["pair_id", "mean", "feature", "endog_dif","mean_rank"]]
    else:
        return None

def getPredValues(dif):
    """
    Calculate the mean of feature means, std for kmers, nullomers, and the endogenous sequence feature means
    
    input
        dif (pd dataframe) - dataframe for one feature, 
            where means of track predictions are summarized for nullomer, control and endogenous sequences

    method
        1. get kmer and nullomer mean of means, 
            as there are multiple nullomer and kmer sequences embedded in the same scaffold,  
            and endogenous mean, 
            for which there is only one sequence and one mean to reflect the track means for that feature
        2. compute standard deviation for kmer, nullomer sequence predictions
        
    return
        (list) - list of summary stat values
    """
    kmer_pred = dif.loc[dif["pair_id"].str.contains("kmer"),"mean"].mean()
    null_pred = dif.loc[dif["pair_id"].str.contains("null"),"mean"].mean()
    ctrl_pred = dif.loc[~dif["pair_id"].str.contains("endog"),"mean"].iloc[0]

    kmer_pred_std = dif.loc[dif["pair_id"].str.contains("kmer"),"mean"].std()
    null_pred_std = dif.loc[dif["pair_id"].str.contains("null"),"mean"].std()
    
    return [kmer_pred, null_pred, ctrl_pred, kmer_pred_std, null_pred_std]
    

def makeLongForm(celldf):
    """
    turn prediction dataframe from df w/ 21907 columns into longform
    
    input
        celldf (pd dataframe) -  sei prediction df w/ 21907 columns
    
    method
        1. use melt to make long form of data, keeping metadata as columns (id_vars) 
        2. split track annotations into cell line, feature, and dataset columns
        
    return 
        test (pd dataframe) - long-form dataframe w/ cell line, feature, and dataset column annotations
        
    """
    # melt dataframe into long form
    test = pd.melt(celldf, id_vars=["nullomer", "coor", "pair_id", "id"],
                           var_name="track", value_name="pred_prob")
    
    # annotate cl, feature, and dataset ids
    test["cl"] = test["track"].apply(lambda x:x.split("|")[0])
    test["feature"] = test["track"].apply(lambda x:x.split("|")[1])
    test["dataset"] = test["track"].apply(lambda x:x.split("|")[2])
    
    return test    
    


def plotHist(dif, coor, feature, ctrl_pred, out):
    """
    plot histogram distribution of mean scores per feature. 
    
    input
        dif (pd dataframe) -  for one feature, prediction means for all nullomer, kmer, and endogenous sequence
        coor (str) - genomic locus
        ctrl_pred (float) - mean feature prediction for endogenous sequence
        out (str) - out file name to use as template for pdf
        
    method 
        1. separate data into types {nullomer, kmer, endog}
        2. for each type, plot histogram distribution
        3. add dashed line for mean of endogenous sequence feature prediction
        4. set title, xlabel, legend
        5. write pdf
    """
    
    #1
    dif["type"] = dif["pair_id"].apply(lambda x: x.split("-")[-1])
    
    #2
    fig, ax = plt.subplots(figsize=(4, 4))

    for t in dif["type"].unique():
        if t!= "endog":
        # plot distribution of prediction differences
            data = dif.loc[dif["type"]==t, "mean"]
            sns.histplot(data, label=t, stat="percent")
    #3
    ax.axvline(ctrl_pred, ls="--", c="g")
    
    #4 name the plot
    ax.set(title=feature,
           xlabel=f"pred - {coor}"
           )
    ax.legend(bbox_to_anchor=(1,1))
    
    #5
    plt.savefig(os.path.join(out.strip(".tsv") + f'.{feature}.pdf'), bbox_inches="tight")
    plt.close()
    

def infoVector(coor, feature, feature_summary_stats):
    """
    make a vector of summary information to write
    
    input
        coor (str) - locus coordinates
        feature (str) - feature name
        feature_summary_stats (list) - list of summary stats. 
        
    method
        1. add meta data to info_vector list
        2. extend list w summary stats
        3. turn every list item into a str
        4. make column list
        
    return
        info (list) - list of meta data and sumamry stats
        columns (list) - names of the data in info list. 
        
    """
    
    # add coordinate and feature names to the list.     
    info_vector = [coor, feature] # prediction track
    
    # add summary stats
    info_vector.extend(feature_summary_stats)
    
    # turn everything into a string
    info = [str(i) for i in info_vector]
    
    # return column names
    columns = ["#coord",
               "track", "predKmer.mean", "predNull.mean",
               "predCtrl", 'predKmer.std',"predNull.std"
               ]

    return info, columns

def main(argv):
    """
    process sei prediction data
    
    filter for:
        - single locus
        - HepG2, K562, and WTC11 cell line feature track predictions
        
    method
        1. Get path variable
            1.1 str split to get coordinates from file name
            1.2 get file row labels
            1.3 make the result .tsv file variables one for cl-specific data, one feature summarized data across all nulls, one for summary of summarized feature data 
        2. if the cl-filtered long form data does not already exist
            2.1 read the hdf5 file
            2.2 turn hdf5 np array data into pd dataframe
            2.3 get the row names
            2.4 format dataframe to add meta data
            2.5 filter for cell line-specific features
            
            2.6 write cell line-specific features longform dataframe
            2.7 else open cldf data
        3. summarize for each feature, turning data into longform
        4. mean of datasets belonging to feature per sequence (nullomer, kmer, endogenous)
        5. compute the means of the nullomer and kmer means
        6. plot histogram if the difference between endogenous and nullomer sequence is greater than 0.05, 
            and the endogenous sequence has predicted probability of at least 0.1 for that feature. 
        7. write summarized means for kmer, nullomer, endogenous sequence per feature. 
        
    """
    

    #1
    PATH = os.path.dirname(FILE)

    #1.1 get locus coordinates from file name
    COOR = (FILE.split("ext4096.")[1]).split("_predictions.h5")[0]
    
    #1.2
    FILE_LABELS = FILE.split("_predictions.h5")[0] + "_row_labels.txt"
    
    #1.3
    OUT_CLDF = os.path.join(PATH, f"HepG2.K562.Pred.{COOR}.tsv")
    OUT_FEATURE_SUMMARY = os.path.join(PATH, f"HepG2.K562.Pred.{COOR}.feature.summary.tsv")
    OUT_SUMMARY = os.path.join(PATH, f"HepG2.K562.Pred.{COOR}.summary.tsv")
    
    #2
    if os.path.exists(OUT_CLDF) is False:
        #2.1 read hdf5 file
        dset = readHdF5(FILE,PATH)

        #2.2 turn data into DF
        df = arrayToDF(dset)

        #2.3 get row names
        names = getRowNames(FILE_LABELS, PATH)

        #2.4 format df
        df = formatDf(df, names, COOR)

        #2.5 keep only data w/ cls
        cldf = filterCL(df)

        #2.6 write cldf
        cldf.drop_duplicates().to_csv(OUT_CLDF, sep='\t', index=False)

        # improve memory efficiency
        del df, dset
        
    else:
        print("already filtered for cldf")
        cldf = pd.read_csv(OUT_CLDF, sep='\t')
    
    #3 compute means, plot hit for locus
    with open(OUT_SUMMARY, "w") as writer:
        
        # make into longform
        longdf = makeLongForm(cldf)
        v = 0
        
        #4 take mean score for every feature dataset, row-wise
        for feature in list(set(longdf["feature"])):

            #4 compute the means, change in means from kmer
            dif = computeFeatureMean(feature, longdf)
            
            # add feature summary data to feature dictionary
            dif.to_csv(OUT_FEATURE_SUMMARY, sep='\t', index=False, header=False, mode="a")
            
            if len(dif)>0:  # some features are empty

                #5 get summary stats about the feature means, ctrl variation
                feature_summary_stats = getPredValues(dif)

                # unpack summary stats variables
                kmer_pred, null_pred, ctrl_pred, kmer_std, null_std = feature_summary_stats

                #6 if nullomer and control are different by atleast 0.05 points, and feature is predicted to have some activity
                if abs(null_pred - ctrl_pred) > 0.05 and ctrl_pred > 0.1:

                    print(feature, COOR, abs(null_pred - ctrl_pred))

                    # plot
                    plotHist(dif, COOR, feature, ctrl_pred, OUT_SUMMARY)

                #7 write data to file
                writeinfo, columns = infoVector(COOR,
                                                feature, feature_summary_stats)
                if v == 0:
                    writer.write("\t".join(columns)+"\n")

                writer.write("\t".join(writeinfo)+"\n")

                v += 1

        writer.close()
        

if __name__ == "__main__":
    main(sys.argv[1:])