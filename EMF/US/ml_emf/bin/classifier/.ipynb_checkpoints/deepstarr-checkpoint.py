import tensorflow as tf

import keras
import keras.layers as kl
from keras.layers.convolutional import Conv1D, MaxPooling1D
from keras.layers.core import Dropout, Reshape, Dense, Activation, Flatten
from keras.layers import BatchNormalization, InputLayer, Input
from keras import models
from keras.models import Sequential, Model
from keras.optimizers import Adam
from keras.callbacks import EarlyStopping, History, ModelCheckpoint

    
import pandas as pd
import numpy as np

from scipy import stats
from scipy.stats import spearmanr
from sklearn.metrics import mean_squared_error, accuracy_score, f1_score

import sys, os
#sys.path.append('Neural_Network_DNA_Demo/')
#from helper # from https://github.com/const-ae/Neural_Network_DNA_Demo
import IOHelper, SequenceHelper 


import random
random.seed(1234)

# SF add
prefix = sys.argv[1] # dataset to pull  
data_path = sys.argv[2]
pred_task = sys.argv[3]
project_dir = os.path.join(data_path, prefix)
if os.path.exists(project_dir) is False:
    os.mkdir(project_dir)
    

model_name="DeepSTARR_ATAC"
params = {'batch_size': 64,
          'epochs': 50, # 100
          'early_stop': 10,
          'kernel_size1': 7,
          'kernel_size2': 3,
          'kernel_size3': 5,
          'kernel_size4': 3,
          'lr': 0.002,
          'num_filters': 256,
          'num_filters2': 60,
          'num_filters3': 60,
          'num_filters4': 120,
          'n_conv_layer': 4,
          'n_add_layer': 2,
          'dropout_prob': 0.4,
          'dense_neurons1': 256,
          'dense_neurons2': 256,
          'pad':'same', 
          "n_nucleotides":4, 
          "seq_len":271,
          "pred_task": pred_task, #"class", # or "reg"
          "n_pred_tasks": 2 # or regression
         }

# SF add params to paramdict
params["prefix"] = prefix
params["data_path"]=data_path
params["project_path"]=project_dir


# SF add change directory
os.chdir(data_path)


# function to load sequences and enhancer activity
def prepare_input(set, prefix):
    # Convert sequences to one-hot encoding matrix
    file_seq = str(f"{prefix}.Sequences_" + set + ".fa")
    input_fasta_data_A = IOHelper.get_fastas_from_file(file_seq, uppercase=True)

    # get length of first sequence
    sequence_length = len(input_fasta_data_A.sequence.iloc[0])

    # Convert sequence to one hot encoding matrix
    seq_matrix_A = SequenceHelper.do_one_hot_encoding(input_fasta_data_A.sequence, sequence_length,
                                                      SequenceHelper.parse_alpha_to_seq)
    print(seq_matrix_A.shape)
    
    X = np.nan_to_num(seq_matrix_A) # Replace NaN with zero and infinity with large finite numbers
    X_reshaped = X.reshape((X.shape[0], X.shape[1], X.shape[2]))

    Activity = pd.read_table(f"{prefix}.Sequences_activity_" + set + ".txt")
    
    Y = []
    for i in Activity.columns[1:]:
        Y.append(Activity[i])
    
    print(set)

    return input_fasta_data_A.sequence, seq_matrix_A, X_reshaped, Y

### Additional metrics

def Spearman(y_true, y_pred):
     return ( tf.py_function(spearmanr, [tf.cast(y_pred, tf.float32), 
                       tf.cast(y_true, tf.float32)], Tout = tf.float32) )

    
def train(selected_model, X_train, Y_train, X_valid, Y_valid, params):

    my_history=selected_model.fit(X_train, Y_train,
                                  validation_data=(X_valid, Y_valid),
                                  batch_size=params['batch_size'], 
                                  epochs=params['epochs'],
                                  callbacks=[EarlyStopping(patience=params['early_stop'], 
                                                           monitor="val_loss", 
                                                           restore_best_weights=True),
                                             History()]
                                 )
    
    return selected_model, my_history
    


# create functions

def summary_statistics(X, Y, set, task):
    pred = main_model.predict(X, batch_size=main_params['batch_size'])
    if task =="US":
        i=0
    if task =="ctrl":
        i=1
        
    if main_params['pred_task'] == "linear":
        print(set + ' MSE ' + task + ' = ' + str("{0:0.2f}".format(mean_squared_error(Y[i], pred[i].squeeze()))))
    
        print(set + ' PCC ' + task + ' = ' + str("{0:0.2f}".format(stats.pearsonr(Y[i], pred[i].squeeze())[0])))
        print(set + ' SCC ' + task + ' = ' + str("{0:0.2f}".format(stats.spearmanr(Y[i], pred[i].squeeze())[0])))
    else:
        
        print(set, "accuracy" + task + " = " + str("{0:0.2f}".format(accuracy_score(Y[i], pred[i].squeeze()))))
        print(set, "f1" + task + " = " + str("{0:0.2f}".format(f1_score(Y[i], pred[i].squeeze()))))

    return pred[i]
    

def DeepSTARR(params=params):
    
    lr = params['lr']
    dropout_prob = params['dropout_prob']
    n_conv_layer = params['n_conv_layer']
    n_add_layer = params['n_add_layer']
    
    # body
    input = kl.Input(shape=(params['seq_len'], params['n_nucleotides']))
    x = kl.Conv1D(params['num_filters'], kernel_size=params['kernel_size1'],
                  padding=params['pad'],
                  name='Conv1D_1st')(input)
    x = BatchNormalization()(x)
    x = Activation('relu')(x)
    x = MaxPooling1D(2)(x)

    for i in range(1, n_conv_layer):
        x = kl.Conv1D(params['num_filters'+str(i+1)],
                      kernel_size=params['kernel_size'+str(i+1)],
                      padding=params['pad'],
                      name=str('Conv1D_'+str(i+1)))(x)
        x = BatchNormalization()(x)
        x = Activation('relu')(x)
        x = MaxPooling1D(2)(x)
    
    x = Flatten()(x)
    
    # dense layers
    for i in range(0, n_add_layer):
        x = kl.Dense(params['dense_neurons'+str(i+1)],
                     name=str('Dense_'+str(i+1)))(x)
        x = BatchNormalization()(x)
        x = Activation('relu')(x)
        x = Dropout(dropout_prob)(x)
    bottleneck = x
    
    # heads per task (developmental and housekeeping enhancer activities)
    
    # SF added below. Accommodate linear and classification tasks.
    pred_task = params["pred_task"] 
    
    if pred_task == "reg":
        activation_ = pred_task
        loss_ = ['mse']*params['n_pred_tasks']
        metrics_ = Spearman
        
    elif pred_task == "class":
        activation_ = tf.nn.softmax
        loss_ = ['binary_crossentropy']*params['n_pred_tasks']
        metrics_ = 'accuracy'
        
    # end SF additions
    
    tasks = ['US', 'ctrl']
    outputs = []
    for task in tasks:
        outputs.append(kl.Dense(1, activation=activation_, name=str('Dense_' + task))(bottleneck))  # changed activation="linear"

    model = keras.models.Model([input], outputs)
    model.compile(keras.optimizers.Adam(lr=lr),
                  loss = loss_,  # SF changed loss=['mse', 'mse'], # loss
                  loss_weights=[1, 1], # loss weigths to balance
                  metrics=[metrics_]) # additional track metric

    return model, params


### MAIN ###

# Data for train/val/test sets
X_train_sequence, X_train_seq_matrix, X_train, Y_train = prepare_input("Train", prefix)
X_valid_sequence, X_valid_seq_matrix, X_valid, Y_valid = prepare_input("Val", prefix)
X_test_sequence, X_test_seq_matrix, X_test, Y_test = prepare_input("Test", prefix)

DeepSTARR()[0].summary()
DeepSTARR()[1] # dictionary

main_model, main_params = DeepSTARR()
main_model, my_history = train(main_model, X_train, Y_train, X_valid, Y_valid, main_params)

# run for each set and enhancer type
train_us_pred = summary_statistics(X_train, Y_train[0], "train", "US")
train_ctrl_pred = summary_statistics(X_train, Y_train[1], "train", "ctrl")
val_us_pred = summary_statistics(X_valid, Y_valid[0], "validation", "US")
val_ctrl_pred = summary_statistics(X_valid, Y_valid[1], "validation", "ctrl")
test_us_pred = summary_statistics(X_test, Y_test[0], "test", "US")
test_ctrl_pred = summary_statistics(X_test, Y_test[1], "test", "ctrl")

preds = pd.merge(pd.DataFrame(test_us_pred), pd.DataFrame(test_ctrl_pred), left_index=True, right_index=True)
preds.columns = ["pred_US", "pred_ctrl"]

                 # write to project dir
os.chdir(project_dir)

# write test predicted and observed values to file
t = pd.merge(pd.DataFrame(Y_test).T, preds, left_index=True, right_index=True)
                 
t.to_csv(f"{model_name}.{prefix}.test.predictions.tsv", sep='\t', index=False)
                 
                 


# write model

## config
model_json = main_model.to_json()
with open(f'Model_{model_name}.{prefix}.json', "w") as json_file:
    json_file.write(model_json)
    
### weights 
main_model.save_weights(f'Model_{model_name}.{prefix}.h5')

# write params # SF added

with open(f'config' + model_name + f'.{prefix}.json', "w") as json_file:
    for key, value in params.items():
        json_file.write(f"{key}:{value}\n")