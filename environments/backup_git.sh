# remove checkpoints
rm ./*/.ipynb_checkpoints

# MPRA flow scripts
echo "backing up mpraflow scripts"
cp -r $HOME/MPRA/mpraflow/bin/ $HOME/wynton_git

# Nullomer crispr
echo "backing up nullomer crispr work"
cp -r $HOME/nullomers/bin-lock/ $HOME/wynton_git/mpra_nullomer_crispr

# for others
echo "backing up for_others" 
cp -r $HOME/other_analyses/for-chengyu_richa/bin/ $HOME/wynton_git/other_analyses/for-chengyu_richa/
cp -r $HOME/other_analyses/for-hai_bats/bin/ $HOME/wynton_git/other_analyses/for-hai_bats/
cp -r $HOME/other_analyses/for-ofer_sei-asthma/bin/ $HOME/wynton_git/other_analyses/for-ofer_sei-asthma/
cp -r $HOME/other_analyses/for-wei_bats/bin/ $HOME/wynton_git/other_analyses/for-wei_bats/

# EMF
echo "backing up EMF" 
cp -r $HOME/EMF/bin/ $HOME/wynton_git/EMF

cp -r $HOME/EMF/US/bin/ $HOME/wynton_git/EMF/US

# bash_profile
echo "backing up bash_profile" 
cp ./.bash_profile $HOME/wynton_git/environments/

# sei
cp -r $HOME/bin/sei-framework/sarah_scripts $HOME/wynton_git/bin/sei-framework/

# mamba env
echo "backing up mamba env" 
conda activate mamba 
conda env export > mamba.yml
mv ./mamba.yml $HOME/wynton_git/environments/
