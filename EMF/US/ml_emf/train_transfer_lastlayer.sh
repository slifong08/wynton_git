#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=40G     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=40G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt=00:20:00   # job requires up to 24 hours of runtime
##$ -t 1-10           # array job with 10 tasks (remove first '#' to enable)
#$ -r y               # if job crashes, it should be restarted
#$ -m ae              # alerts to mail about
#$ -M sarah.fong@ucsf.edu # email to mail to

#$ -e /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".error.txt
#$ -o /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".train_transfer.output.txt

## If you array jobs (option -t), this script will run T times, once per task.
## For each run, $SGE_TASK_ID is set to the corresponding task index (here 1-10).
## To configure different parameters for each task index, one can use a Bash 
## array to map from the task index to a parameter string.

## All possible parameters
# params=(1bac 2xyz 3ijk 4abc 5def 6ghi 7jkl 8mno 9pqr 10stu)

## Select the parameter for the current task index
## Arrays are indexed from 0, so we subtract one from the task index
# param="${params[$((SGE_TASK_ID - 1))]}"

date
hostname
ml load cuda/11.5 gcc/7.3.1

source activate /wynton/home/ahituv/fongsl/.conda/envs/torch

#python3 /wynton/home/ahituv/fongsl/EMF/US/ml_emf/bin/train_transfer_lastlayer.py --train_valid_path "$1" --foldify --delimiter tab --seed 42 --train_batch_size 1024 --train_workers 8 --valid_batch_size 4098 --valid_workers 8 --epoch_num 10 --batch_per_epoch 500 --weights uniform --seqsize "$2" --temp .TEMPDIR --use_single_channel --singleton_definition integer --gpu 0 --model_dir "$3" --ks 7 --blocks 256 128 128 64 64 64 64 --resize_factor 4 --se_reduction 4 --shift 0.5 --scale 0.5 --loss kl --final_ch 18 --optimizer adamw --model /wynton/home/ahituv/fongsl/EMF/US/ml_emf/data/legnet/model_300.pth "$4" "$5"

# turn arguments into a str
str="$*"

# echo
echo python3 /wynton/home/ahituv/fongsl/EMF/US/ml_emf/bin/train_transfer_lastlayer.py $str

# run
python3 /wynton/home/ahituv/fongsl/EMF/US/ml_emf/bin/train_transfer_lastlayer.py $str

# where
# $1 is training/validation data
# $2 is sequence size
# $3 is the model dir

# need variables - number of epochs? batch per epoch size? 

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"
