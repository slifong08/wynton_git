#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=1G     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=2G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt=24:00:00   # job requires up to 24 hours of runtime

#$ -r y               # if job crashes, it should be restarted
#$ -m ae              # alerts to mail about
#$ -M sarah.fong@ucsf.edu # email to mail to
#$ -e /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".error.txt
#$ -o /wynton/home/ahituv/fongsl/qsub/log/"$JOB_ID".output.txt


echo date
echo hostname

#micromamba activate /wynton/home/ahituv/fongsl/.micromamba/envs/mamba # environment
cd /wynton/group/ahituv/bin
       

"""# index
java -Xmx6000m -jar FlashFry-assembly-1.15.jar index \
    --tmpLocation /wynton/group/ahituv/fongsl/projects/nullomers/data/flashfry/tmp \
    --database hs1_database \
    --reference /wynton/group/ahituv/data/dna/hs1/hs1.fa \
    --enzyme spcas9ngg
    
"""
# discover
java -Xmx4g -jar FlashFry-assembly-1.15.jar \
discover \
 --database ./hs1_database \
 --fasta /wynton/group/ahituv/fongsl/projects/CRISPR_nullomer_MPRA/data/20230815_nullomer_MPRA_assoc/15mer.fo.pam.scaffold.ext200.library.TWIST.noAdaptors.fa \
 --output /wynton/group/ahituv/fongsl/projects/CRISPR_nullomer_MPRA/data/flashfry/15mer.fo.pam.scaffold.ext200.library.TWIST.noAdaptors.fa.sites
 

# score
java -Xmx4g -jar FlashFry-assembly-1.15.jar \
score \
 --input /wynton/group/ahituv/fongsl/projects/CRISPR_nullomer_MPRA/data/flashfry/15mer.fo.pam.scaffold.ext200.library.TWIST.noAdaptors.fa.sites \
 --output /wynton/group/ahituv/fongsl/projects/CRISPR_nullomer_MPRA/data/flashfry/15mer.fo.pam.scaffold.ext200.library.TWIST.noAdaptors.fa.sites.scored \
 --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot,bedannotator \
 --database ./hs1_database \
--inputAnnotationBed /wynton/group/ahituv/fongsl/projects/CRISPR_nullomer_MPRA/data/20230815_nullomer_MPRA_assoc/15mer.fo.pam.scaffold.ext200.library.TWISThs1.bed \
--transformPositions /wynton/group/ahituv/fongsl/projects/CRISPR_nullomer_MPRA/data/20230815_nullomer_MPRA_assoc/15mer.fo.pam.scaffold.ext200.library.TWISThs1.bed



## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"