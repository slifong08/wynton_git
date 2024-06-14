#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=12G     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=12G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt=00:02:00   # job requires up to 24 hours of runtime

#$ -r y               # if job crashes, it should be restarted
#$ -m ae              # alerts to mail about
#$ -M sarah.fong@ucsf.edu # email to mail to
#$ -e /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".error.txt
#$ -o /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".output.txt


echo date
echo hostname
ml load openjdk/17

# add miniconda from bin
#echo init miniconda3 from /wynton/group/ahituv/bin/
#/wynton/home/ahituv/fongsl/micromamba/bin/micromamba init bash

# source the miniconda from local directory
#echo sourcing miniconda3 from $home/.bashrc
#source $HOME/.bashrc


# activate virtual environment
#echo activate virtual environment
micromamba activate $HOME/micromamba/envs/env_nf #

# add symlink because nextflow is a piece of work
#echo add symlink for miniconda activate
#ln -s /wynton/group/ahituv/bin/miniconda3/bin/activate $HOME/.conda/envs/MPRAflow/bin/activate


# export conda cache dir
#export NXF_CONDA_CACHEDIR=/wynton/group/ahituv/fongsl/src/MPRAflow/work/conda

# go to directory
echo change dir /wynton/group/ahituv/bin/pipelines/nf-core-atacseq-dev/workflow/
cd /wynton/group/ahituv/bin/pipelines/nf-core-atacseq-dev/workflow/


###
# NEXTFLOW COMMAND
###
#micromamba activate /wynton/home/ahituv/fongsl/micromamba/envs/atac 
#micromamba


#echo nextflow run nf-core/atacseq --input /wynton/group/ahituv/fongsl/projects/EMF/data/ATAC/WTC11NGN2.1/samples.csv --outdir /wynton/group/ahituv/fongsl/projects/EMF/data/ATAC/WTC11NGN2.1/results --genome GRCh38 --read_length 150  -profile singularity

nextflow run /wynton/group/ahituv/bin/pipelines/nf-core-atacseq-dev/workflow/ --input /wynton/group/ahituv/fongsl/projects/EMF/data/ATAC/WTC11NGN2.1/samples.csv --outdir /wynton/group/ahituv/fongsl/projects/EMF/data/ATAC/WTC11NGN2.1/results --fasta /wynton/group/ahituv/data/dna/hg38/hg38.fa.gz --gtf /wynton/group/ahituv/data/dna/hg38/hg38.knownGene.gtf.gz --genome GRCh38 --read_length 150  -profile singularity --narrow_peak -resume

#nextflow run nf-core/atacseq --input /wynton/group/ahituv/fongsl/projects/EMF/data/ATAC/WTC11NGN2.1/samples.csv --outdir /wynton/group/ahituv/fongsl/projects/EMF/data/ATAC/WTC11NGN2.1/results --genome GRCh38 --read_length 150  -profile singularity 

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"