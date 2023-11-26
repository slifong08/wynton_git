#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l h_vmem=30G      # memory 
#$ -l mem_free=30G     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=30G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt=36:15:00   # job requires up to 24 hours of runtime
##$ -t 1-10           # array job with 10 tasks (remove first '#' to enable)
#$ -r y               # if job crashes, it should be restarted
#$ -m ae              # alerts to mail about
#$ -M sarah.fong@ucsf.edu # email to mail to
#$ -e /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".error.txt
#$ -o /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".output.txt
##$ -t 1-2:1        # job range:step size
##$ -tc 1             # n jobs to run at once

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

# ml load
ml load Sali miniconda3/23.3.1-0-py39 #CBI picard/2.27.5  samtools/1.18  bwa/0.7.17

# activate env
source activate /wynton/home/ahituv/fongsl/.conda/envs/mpraflow_wynton # environment

# change dir to datapath
cd "$1"

# step1 - make config
echo /wynton/home/ahituv/fongsl/MPRA/mpraflow/bin/make_mpraflow_config.py
python3 /wynton/home/ahituv/fongsl/MPRA/mpraflow/bin/make_mpraflow_config.py "$1" "$2" "$3" "$4" "$5" "$6" 

echo  /wynton/home/ahituv/fongsl/MPRA/mpraflow/bin/association.py
python3  /wynton/home/ahituv/fongsl/MPRA/mpraflow/bin/association.py "$1" "$2"

# $1 is .bed file to intersect

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"

"""
Args for make_mpraflow_config.py

$1 arg_parser.add_argument("path", type=str, help='full path to fastq read directory')
$2 arg_parser.add_argument("name", type=str, help='name of experiment')
$3 arg_parser.add_argument("INS", type=str, help='Fastq insert reads')
$4 arg_parser.add_argument("PE",  type=str, help='Fastq Paired end reads')
$5 arg_parser.add_argument("barcode", type=str, help='barcodes as demultiplexed i7 fastq reads')
$6 arg_parser.add_argument("design", type=str, help='fasta file of library design')


# set directory w/ scripts

BIN = "HOME"
echo $BIN
if [ $BIN = "HOME" ]
then
    BINDIR = "/wynton/home/ahituv/fongsl/MPRA/mpraflow/bin"
else
    BINDIR = "/wynton/group/ahituv/fongsl/src/MPRAflow/association/"
fi
"""
"""