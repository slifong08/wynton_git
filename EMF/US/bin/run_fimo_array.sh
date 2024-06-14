#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=2G     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=2G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt=00:10:00   # job requires up to 24 hours of runtime
##$ -t 1-1876:1        # job range:step size
##$ -tc 20            # n jobs to run at once
#$ -r y               # if job crashes, it should be restarted
#$ -m ae              # alerts to mail about
#$ -M sarah.fong@ucsf.edu # email to mail to

#$ -e /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".error.txt
#$ -o /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".train_transfer.output.txt



date
hostname
#ml load cuda/11.5 #gcc/7.3.1

micromamba activate /wynton/home/ahituv/fongsl/micromamba/envs/mamba

#cd 

# run fimo on ultrasound library

python3 $HOME/EMF/US/bin/run_fimo.py "$1" "$2" $SGE_TASK_ID
# 1 outdir
# 2 fasta

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"
