# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# load anaconda
if [[ -n "$MODULEPATH" ]]
then
    module load Sali
    #module load CBI miniconda3-py39/4.12.0
    module load CBI miniconda3/23.3.1-0-py39

fi
# User specific environment and startup programs

PATH=$PATH:$HOME/.local/bin:$HOME/bin
export PATH

# python sarah custom tools path
PYTHONPATH=$PYTHONPATH:$PWD/tools/py_
export PYTHONPATH

#python sarah custom tools genome
PYTHONPATH=$PYTHONPATH:$HOME/tools/genome
export PYTHONPATH

# pytho sarah custom tools evo 
PYTHONPATH=$PYTHONPATH:$HOME/tools/evo
export PYTHONPATH

# JAVA v 17
JAVA_HOME=/wynton/home/ahituv/fongsl/bin/lib/jvm/jdk-17.0.6
export JAVA_HOME

# UCSC tools
export PATH=$PATH:$HOME/bin/ucsc_exe

# memesuite
export PATH=/wynton/home/ahituv/fongsl/bin/meme-5.5.1/bin:/wynton/home/ahituv/fongsl/meme/libexec/meme-5.5.1:$PATH

# argtable
export PATH=$PATH:$HOME/bin/argtable2

# github cli
export GH_HOST='slifong08'
export GH_ENTERPRISE_TOKEN='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'

# tunneling from wynton - alias
alias 'nb1'="jupyter notebook --no-browser --ip='*' --port=7778"
alias 'nb2'="jupyter notebook --no-browser --ip='*' --port=7779"

# tunnel from wynton - command
function tnl { jupyter notebook --no-browser --ip='*' --port=$1; }
export -f tnl

# set up a new project folder
function project_setup { python ./tools/py_/setup_new_project.py $1; }
export -f project_setup

# environments
alias env="source activate mamba"
alias sei="source activate sei"

# wynton short-cuts
alias fsl="cd $HOME"
alias dev2="ssh dev2.wynton.ucsf.edu"
alias temp='cd $TMPDIR'

# shared resources dirs
alias group="cd /wynton/group/ahituv"
alias groupfsl="cd /wynton/group/ahituv/fongsl"
alias dc="cd /wynton/group/databases/"

# trash commands
alias rm="/wynton/home/ahituv/fongsl/tools/conda_env/trash.pl"
alias rrm="'rm'"
alias trash="rrm -r /wynton/scratch/fongsl/Trash/"

# project paths
alias lock="cd $HOME/nullomers/data/lock"
alias dna="cd $HOME/dna/"
alias us='cd $HOME/EMF/US'

# box downloads
alias box='lftp --user sarah.fong@ucsf.edu ftps://ftp.box.com'

# modules
alias cuda='module load cuda/11.5'

# sequence alignment toolkit
alias sequencing='ml load bcftools/1.18 bcl2fastq/2.20.0 bowtie/1.3.1 bowtie2/2.5.1 bwa/0.7.17 hisat2/2.2.0'

# ncbi toolkit
alias ncbi='ml load sratoolkit/3.0.0 blast/2.14.0 blast+/2.12.0 blat/37x1'
