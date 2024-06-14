#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l h_vmem=30G     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=20G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt=02:30:00   # job requires up to 24 hours of runtime
##$ -t 1-10           # array job with 10 tasks (remove first '#' to enable)
#$ -r y               # if job crashes, it should be restarted
#$ -m ae              # alerts to mail about
#$ -M sarah.fong@ucsf.edu # email to mail to

#$ -e /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".error.txt
#$ -o /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".deepstarr.cpu.output.txt

date
hostname
#ml load cuda/11.5 #gcc/7.3.1

micromamba activate /wynton/home/ahituv/fongsl/micromamba/envs/DeepSTARR

# go to home
cd /wynton/home/ahituv/fongsl/EMF/US/ml_emf/bin/classifier

#export CUDA_VISIBLE_DEVICES=$SGE_GPU

#echo sgegpu: $SGE_GPU
echo training dataset: $1
echo dir: $2
echo prediction task: $3
echo standard scale outputs: $4


python ./deepstarr.py $1 $2 $3 $4

# 1 = prefix
# 2 = data_path
# 3 =  prediction task "class" | "reg"
# 4 = standard scaling? Bool.
# 5 = n prediction tasks

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"
