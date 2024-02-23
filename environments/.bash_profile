# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs

#ml load CBI emacs/29.1 miniconda3/23.5.2-0-py311


# gpu
alias cuda="ml load cuda/11.5"
alias torch="conda activate torch"

# working virtual env
alias env="conda activate mamba"


# tunneling from wynton - alias
alias 'nb1'="python3 -m notebook --no-browser --ip='*' --port=7778"
alias 'nb2'="python3 -m notebook --no-browser --ip='*' --port=7779"

# tunnel from wynton - command
function tnl { python3 -m notebook --no-browser --ip='*' --port=$1; }
export -f tnl

# set up a new project folder
function project_setup { python $HOME/tools/py_/setup_new_project.py $1; }
export -f project_setup

# install new ipython kernels
function kernel_install { python3 -m ipykernel install --user --name $1 --display-name $1; }
export -f kernel_install

# environments
alias env="source activate mamba"
alias sei="source activate sei"

# wynton short-cuts
alias fsl="cd $HOME"
alias dev2="ssh dev2.wynton.ucsf.edu"
alias temp='cd $TMPDIR'
alias scratch='cd /scratch/fongsl/'
alias mamba='micromamba'

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
alias usml='cd $HOME/EMF/US/ml_emf/bin/human_legnet'
# box downloads
alias box='lftp --user sarah.fong@ucsf.edu ftps://ftp.box.com'

# bash/python conventions
alias python="python3"
alias run="."

# python paths
# sarah custom tools path
PYTHONPATH=$PYTHONPATH:$PWD/tools/py_
export PYTHONPATH

# sarah custom tools genome
PYTHONPATH=$PYTHONPATH:$HOME/tools/genome
export PYTHONPATH

# sarah custom tools evo 
PYTHONPATH=$PYTHONPATH:$HOME/tools/evo
export PYTHONPATH

# export local path
export PATH=$PATH:/wynton/home/ahituv/fongsl/.local/bin
