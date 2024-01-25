# -*- coding: utf-8 -*-
from argparse import ArgumentParser
from collections import Counter
import json
import numpy as np
import pandas as pd
from pathlib import Path
import random 
from scipy.stats import pearsonr, spearmanr
import scipy

import torch
import torch.nn as nn
from torch.utils.data import Dataset
import torch.nn.functional as F 

from torch.utils.data import DataLoader
from torch.utils.data import Sampler
from torch import Tensor 
from typing import Sequence, Iterator, ClassVar

from model import SeqNN, TLModel

torch.cuda.empty_cache()
###
# parse arguments
###
parser = ArgumentParser()
TRAIN_VAL_PATH = "/wynton/home/ahituv/fongsl/EMF/US/ml_emf/data/legnet/complex_Native.txt"

parser.add_argument("--valid_folds", nargs='+',  type=int, default=list())
parser.add_argument("--seed", type=int, default=42)
parser.add_argument("--valid_batch_size", type=int, default=4098)
parser.add_argument("--valid_workers", type=int, default=8)
parser.add_argument("--batch_per_epoch", type=int, default=1000)
parser.add_argument("--seqsize", type=int, default=120)
parser.add_argument("--temp", default=".TEMPDIR", type=Path)
parser.add_argument("--use_single_channel", action="store_true")
parser.add_argument("--singleton_definition", choices=["integer", "threshold1100"], default="integer")
parser.add_argument("--gpu", type=int, default=0)
parser.add_argument("--ks", default=5, type=int, help="kernel size of convolutional layers")
parser.add_argument("--blocks", default=[256, 256, 128, 128, 64, 64, 32, 32], nargs="+", type=int, help="number of channels for EffNet-like blocks")
parser.add_argument("--resize_factor", default=4, type=int, help="number of channels in a middle/high-dimensional convolutional layer of an EffNet-like block")
parser.add_argument("--se_reduction", default=4, type=float, help="reduction number used in SELayer")
parser.add_argument("--loss", choices=["mse", "kl"], default="mse", type=str, help="loss function")
parser.add_argument("--final_ch", default=18, type=int, help="number of channels of the final convolutional layer")
parser.add_argument("--target", type=str, help="path to a testing dataset")
parser.add_argument("--model", type=str, default="/wynton/home/ahituv/fongsl/EMF/US/ml_emf/data/legnet/model_80.pth", help="path to a .pth file where parameters are stored")
parser.add_argument("--output", type=str, default="results.tsv", help="path to file where results will be stored")
parser.add_argument("--output_format", type=str, choices=["tsv", "csv", "json"], default="tsv", help="foramt of the output file")
parser.add_argument("--delimiter", default='space', type=str, help="delimiter that separates columns in a training file")
parser.add_argument("--tl_layers", action='store_true', help="add new input/output layers")

args = parser.parse_args()



def preprocess_data(data, length):
    # deal with plasmid # 
    # changed so that there is no scaffold to insert test sequence into, just 500 N's
    # changed so that there is no adaptor adjustment. Assumes sequences DO NOT have adaptor

    LEFT_ADAPTER, RIGHT_ADAPTER = "", ""

    PLASMID = "N" * 500
    PLASMID = PLASMID.upper()

    INSERT_START = PLASMID.find('N'*80)
    
    print("PLASMID", PLASMID, "INSERT_START", INSERT_START)
    data = data.copy()
    
    add_part = PLASMID[INSERT_START-length:INSERT_START]  # find where the sequence insert goes, moot point in N plasmid
    print("ADD_PART", add_part, len(add_part))
    
    data.seq = data.seq.apply(lambda x:  add_part + x[len(LEFT_ADAPTER):]) # add sequence to scaffold, excluding left adaptor    
    print("after lambda:", len(data["seq"].iloc[0]), data["seq"].iloc[0])

    data.seq = data.seq.str.slice(-length, None) # remove the right side of the sequence. 
    print("after splice:", len(data["seq"].iloc[0]), data["seq"].iloc[0])
    return data  # return dataframe w/ sequences of 150 length

def set_global_seed(seed: int) -> None:
    """
    Sets random seed into PyTorch, TensorFlow, Numpy and Random.
    Args:
        seed: random seed
    """
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.backends.cudnn.deterministic = True #type: ignore
    torch.backends.cudnn.benchmark = False #type: ignore


def n2id(n):
    """one-hot-like dict, seq -> code"""
    
    CODES = {
        "A": 0,
        "T": 3,
        "G": 1,
        "C": 2,
        'N': 4
    }

    return CODES[n.upper()]

def id2n(i):
    """reverse one-hot-like dict, code -> seq"""
    CODES = {
        "A": 0,
        "T": 3,
        "G": 1,
        "C": 2,
        'N': 4
    }

    INV_CODES = {value: key for key, value in CODES.items()}
    return INV_CODES[i]

def n2compl(n):
    
    """ seq -> complementary seq"""
    
    COMPL = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G',
    'N': 'N'
    }
    
    return COMPL[n.upper()]

class Seq2Tensor(nn.Module):
    '''
    Encode sequences using one-hot encoding after preprocessing.
    '''
    def __init__(self):
        super().__init__()

    def forward(self, seq: str) -> torch.Tensor:
        seq_i = [n2id(x) for x in seq]
        code = torch.from_numpy(np.array(seq_i))
        code = F.one_hot(code, num_classes=5) # 5th class is N, one_hot encodeing function from pytorch

        code = code[:, :5].float()  # turn everything into a float value
        code[code[:, 4] == 1] = 0.25 # encode Ns with .25
        code =  code[:, :4]
        return code.transpose(0, 1)  # return transpose of the code

def parameter_count(model):
    pars = 0  
    for _, p  in model.named_parameters():    
        pars += torch.prod(torch.tensor(p.shape))
    return pars

class SeqDatasetRev(Dataset):
    def __init__(self, ds_rev, size, use_single_channel):
        self.ds = ds_rev
        self.size = size
        self.use_single_channel = use_single_channel
        self.totensor = Seq2Tensor() 
        
    def transform(self, x):
        assert isinstance(x, str)
        print("ASSERTION ERROR?", len(x), self.size)
        assert len(x) == self.size
        return self.totensor(x)
    
    def __getitem__(self, i):
        seq = self.transform(self.ds.seq.values[i])  # check that sequence is specified size
        rev = torch.full( (1, self.size), self.ds.rev.values[i], dtype=torch.float32) # load tensor with sequence values
        
        if self.use_single_channel:
            single = torch.full( (1, self.size) , self.ds.is_singleton.values[i], dtype=torch.float32)
            X = torch.concat([seq, rev, single])  # combine
        else:
            X = torch.concat([seq, rev], dim=0)
        
        bin = self.ds.bin.values[i]
        return X, bin 
    
    def __len__(self):
        return len(self.ds.seq)
    
POINTS = np.array([-np.inf, *range(1, 18, 1), np.inf])

class SeqDatasetRevProb(Dataset):

    """Assign sequence class probability"""
    
    def __init__(self, ds_rev, size, use_single_channel, shift=0.5, scale=1.5):
        self.ds = ds_rev
        self.size = size
        self.totensor = Seq2Tensor() 
        try:
            self.scale = float(scale)
        except ValueError:
            self.scale = scale
            print("Using adaptive scale")
        self.shift = shift 
        self.use_single_channel = use_single_channel
        
    def transform(self, x):
        assert isinstance(x, str)
        #print("ASSERTION ERROR?", len(x), self.size)
        assert len(x) == self.size
        return self.totensor(x)
    
    def __getitem__(self, i):
        seq = self.transform(self.ds.seq.values[i])
        rev = torch.full( (1, self.size), self.ds.rev.values[i], dtype=torch.float32)
        if self.use_single_channel:
            single = torch.full( (1, self.size) , self.ds.is_singleton.values[i], dtype=torch.float32)
            X = torch.concat([seq, rev, single], dim=0)
        else:
            X = torch.concat([seq, rev], dim=0)
        bin = self.ds.bin.values[i]
        if isinstance(self.scale, float):
            norm = scipy.stats.norm(loc=bin + self.shift, scale=self.scale)
        elif self.scale == "adaptive":
            s = -0.1364 * bin + 2.7727
            norm = scipy.stats.norm(loc=bin + self.shift, scale=s)
        else:
            raise Exception("Wrong scale")
        
        cumprobs = norm.cdf(POINTS)
        probs = cumprobs[1:] - cumprobs[:-1]
        return X, probs, bin
    
    def __len__(self):
        return len(self.ds.seq)

class CustomWeightedRandomSampler(Sampler[int]):
    r"""Samples elements from ``[0,..,len(weights)-1]`` with given probabilities (weights).
    Args:
        weights (sequence): a sequence of weights, not necessary summing up to one
        num_samples (int): number of samples to draw
        replacement (bool): if ``True``, samples are drawn with replacement.
            If not, they are drawn without replacement, which means that when a
            sample index is drawn for a row, it cannot be drawn again for that row.
        generator (Generator): Generator used in sampling.
    Example:
        >>> list(WeightedRandomSampler([0.1, 0.9, 0.4, 0.7, 3.0, 0.6], 5, replacement=True))
        [4, 4, 1, 4, 5]
        >>> list(WeightedRandomSampler([0.9, 0.4, 0.05, 0.2, 0.3, 0.1], 5, replacement=False))
        [0, 1, 4, 3, 2]
    """
    weights: Tensor
    num_samples: int
    replacement: bool

    SAMPLES_PER_GROUP: ClassVar[int] = 2 ** 12

    def __init__(self, weights: Sequence[float], num_samples: int,
                 replacement: bool = True, generator=None) -> None:
        if not isinstance(num_samples, int) or isinstance(num_samples, bool) or \
                num_samples <= 0:
            raise ValueError("num_samples should be a positive integer "
                             "value, but got num_samples={}".format(num_samples))
        if not isinstance(replacement, bool):
            raise ValueError("replacement should be a boolean value, but got "
                             "replacement={}".format(replacement))    
        
        self.weights = torch.as_tensor(weights, dtype=torch.double)
        self.num_samples = num_samples
        self.replacement = replacement
        self.generator = generator
        self.groups, self.group_weights, self.inner_weights = self.groups_split()

    def __iter__(self) -> Iterator[int]:
        group_ids = torch.multinomial(self.group_weights, self.num_samples, replacement=True, generator=self.generator)
        gr, counts = torch.unique(group_ids , sorted=True, return_counts=True)
        smpls = []
        for i in range(len(gr)):
            g = gr[i]
            rand_tensor = torch.multinomial(self.inner_weights[g], counts[i], self.replacement, generator=self.generator)
            idxs = self.groups[g][0] + rand_tensor
            smpls.append(idxs)
        rand_tensor = torch.concat(smpls)
        pi = torch.randperm(rand_tensor.shape[0])
        rand_tensor = rand_tensor[pi] 
        yield from iter(rand_tensor.tolist())

    def __len__(self) -> int:
        return self.num_samples

    def groups_split(self):
        groups = []
        N = len(self.weights)
        for start in range(0, N, self.SAMPLES_PER_GROUP):
            end = min(start+self.SAMPLES_PER_GROUP, N)
            groups.append(torch.arange(start, end, 1, dtype=torch.long))
        inner_weights = [self.weights[ids] for ids in groups]
        weights = torch.FloatTensor([w.sum() for w in inner_weights])
        return groups, weights, inner_weights

def get_weights(df, tp):
    if tp == "uniform":
        weights = np.full_like(df.bin.values, fill_value=1)
        return weights
    elif tp == "counts":
        weights = df.cnt.values
        return weights
    raise NotImplementedError()


class DataloaderWrapper:
    def __init__(self, dataloader, batch_per_epoch):
        self.batch_per_epoch = batch_per_epoch
        self.dataloader = dataloader
        self.iterator = iter(dataloader)

    def __len__(self):
        return self.batch_per_epoch
    
    def __next__(self):
        try:
            return next(self.iterator)
        except StopIteration:
            self.iterator = iter(self.dataloader)

    def __iter__(self):
        for _ in range(self.batch_per_epoch):
            try:
                yield next(self.iterator)
            except StopIteration:
                self.iterator = iter(self.dataloader)

def revcomp(seq):
    """get reverse complement sequence w n2compl"""
    return "".join((n2compl(x) for x in reversed(seq)))


def get_rev(df):
    """copy, transform dataframe w/ reverse complement sequence"""
    revdf = df.copy()
    revdf['seq'] = df.seq.apply(revcomp)
    revdf['rev'] = 1
    return revdf
    

def add_rev(df):
    """concat forward and reverse complement dataframes"""
    df = df.copy()
    revdf = df.copy()
    revdf['seq'] = df.seq.apply(revcomp)
    df['rev'] = 0
    revdf['rev'] = 1
    df = pd.concat([df, revdf]).reset_index(drop=True)
    return df


def infer_singleton(arr, method):
    """Method about singleton inference. Still not sure what a singleton is"""
    
    if method == "integer":
        return np.array([x.is_integer() for x in arr])  # make a np array of integer values
    
    elif method.startswith("threshold"):
        th = float(method.replace("threshold", ""))  # make flow array for instances above a threshold value.
        cnt = Counter(arr)
        return np.array([cnt[x] >= th for x in arr])
    else:
        raise Exception("Wrong method")


def add_singleton_column(df, method):
    """determine whether there is only one value in bin? i.e. is singleton."""
    
    df = df.copy()
    df["is_singleton"] = infer_singleton(df.bin.values,method)
    return df 
            

### 
# Load model
###

print("Loading model...")
device = torch.device(f"cuda:{args.gpu}")  # get device (CUDA or CPU) 
print("cuda?", torch.cuda.is_available())
print(device)

###
# instantiate model
###

model = SeqNN(seqsize=args.seqsize, use_single_channel=args.use_single_channel, block_sizes= args.blocks, ks=args.ks, 
              resize_factor=args.resize_factor, se_reduction=args.se_reduction, final_ch=args.final_ch).to(device)

if args.tl_layers is True:
    # transfer learning model architecture
    model = TLModel(model, seqsize=args.seqsize, use_single_channel=args.use_single_channel, block_sizes=args.blocks, ks=args.ks,
              resize_factor=args.resize_factor, se_reduction=args.se_reduction, final_ch=args.final_ch).to(device)  

# load model
model.load_state_dict(torch.load(args.model, map_location=device))
model.eval()


# measure params
print('Parameters:', int(parameter_count(model)))

###
# load sequence data
###

# load input sequences
print(f"Reading dataset {args.target}...")
folds = args.valid_folds
if folds:
    cols = ['seq', 'bin']
else:
    cols = ['seq', 'bin', 'fold']

# load data
df = pd.read_table(args.target, 
     sep='\t' if args.delimiter == 'tab' else ' ', 
     header=None) # modified

# Does sequence data have ground truth? 
if len(df.columns) == 1:
    df["bin"] = 1.0  # 20231016 SF: added as dummy column if no bin column
    print(f'Warning: there is just one column in the dataframe.\
            Are you sure that you provided the correct delimiter?\
            Using {args.delimiter} now.')
    
# rename columns depending on number of columns in dataframe. 
df.columns = ['seq', 'bin', 'fold'][:len(df.columns)]
print(df.head())

# preprocess data by inserting 80mer into 150bp promoter scaffold
df_orig = df.copy()
df = preprocess_data(df, args.seqsize)

# use single_channel?
if args.use_single_channel:
    df = add_singleton_column(df, args.singleton_definition)

# split data into folds
if folds:
    df = df[df.fold.isin(folds)]
    
#df = add_rev(df)

# get reverse of sequence for data augmentation
df_rev = get_rev(df)
df['rev'] = 0  # binary column signaling that this sequence is not the reverse sequence. 

#tmp_filename = args.target + '.tmp'

# make tensors for forward and reverse
ds = SeqDatasetRev(df, size=args.seqsize, use_single_channel=args.use_single_channel)
ds_rev = SeqDatasetRev(df_rev, size=args.seqsize, use_single_channel=args.use_single_channel)

# load data? Not sure where this function goes to. 
dl = DataLoader(ds, 
                batch_size=args.valid_batch_size, 
                num_workers=args.valid_workers,
                shuffle=False)

dl_rev = DataLoader(ds_rev, 
                    batch_size=args.valid_batch_size, 
                    num_workers=args.valid_workers,
                    shuffle=False)

# make predictions!
print("Evaluation...")
y, y_true, y_rev = list(), list(), list()

# for forward sequence
for its in dl:
    #print("ITS", its)
    if type(its) in (tuple, list):
        x, yt = its
        x = x.to(device)
        y.extend(model(x)[-1].detach().cpu().flatten().tolist()) # here, predict sequence vector in model 
        y_true.extend(yt.detach().cpu().flatten().tolist())  # add y_true value to vector
    else:
        its = its.to(device)
        y.extend(model(its)[-1].detach().cpu().flatten().tolist())  # just make predictions

# for reverse sequences
for its in dl_rev:
    #print("REV_ITS", its)
    if type(its) in (tuple, list):
        x, yt = its
        x = x.to(device)
        y_rev.extend(model(x)[-1].detach().cpu().flatten().tolist())
    else:
        its = its.to(device)
        y_rev.extend(model(its)[-1].detach().cpu().flatten().tolist())
        
# check that predictions are the same
assert len(y) == len(y_rev)


# take mean prediction of forward and reverse? 
for i in range(len(y)):
    y[i] = (y[i] + y_rev[i]) / 2
    
# write prediction output
if args.output_format == 'json':
    keys = list(range(0, len(y)))
    d = dict(zip(keys, y))
    with open(args.output, 'w') as f:
        json.dump(d, f)
else:
    df = df_orig[['seq', 'bin']]
    df['bin'] = y # assign the predictions 
    df.to_csv(args.output, sep='\t' if args.output_format == 'tsv' else ',', index=None, header=False)

# measure performance of predictions    
y = np.array(y)
if y_true:
    try:
        y_true = np.array(y_true)
        mse = np.mean((y_true - y) ** 2)
        r_pearson = pearsonr(y, y_true)
        r_spearman = spearmanr(y, y_true)
        print(f'MSE: {mse}, Pearson: {r_pearson}, Spearman: {r_spearman}')
    except ValueError:
        pass