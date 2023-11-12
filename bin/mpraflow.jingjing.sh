1.Association #### Change NUM, $$BC$$,ADAPTOR  accordingly ## where
"""
ADAPTOR (str) is the sequence adaptors? -- yes
1:53
NUM1 = length(R1+BC)+1 (edited) 
1:54
NUM2 = length(BC)
1:55
NUM3 = length(R2) (edited) 
BC - string for 1 barcode (e.g. AGTGCGTAGCTAGTC)
"""

paste <(zcat ../../barcode/${SAMPLE}*R1_001.fastq.gz) <(zcat ../../barcode/${SAMPLE}*R3_001.fastq.gz) <(zcat ../../barcode/${SAMPLE}*R2_001.fastq.gz) | awk '{ count+=1; if ((count == 1) || (count == 3)) { print $1 } else { print $1"$$BC$$"$2$3 }; if (count == 4) { count=0 } }' > ${SAMPLE}.fastq 

python2.7 SplitFastQdoubleIndexBAM.py -s NUM1 -l NUM2 -m NUM3 -i ../index.lst --remove --summary ../${SAMPLE}.fastq
python2.7 MergeTrimReadsBAM.py --mergeoverlap -f FORWARD-ADAPTOR -s REVERSE-ADAPTOR demultiplex.bam 

samtools view -F -r ${SAMPLE}/MergeTrim.bam | awk 'BEGIN{ OFS= "\t" }{ for (i=12; i<=NF; i++) { if ($i ~ /^XJ:Z:/) print $10,substr($i,6) }}' | sort | gzip -c  > raw_counts/${SAMPLE}.tsv.gz
zcat raw_counts/${SAMPLE}.tsv.gz  | grep -v "N" | awk 'BEGIN{ OFS="\t" }{ if (length($1) == 15) { print } }' | gzip -c > filter_counts/${SAMPLE}.filtered.tsv.gz

bwa mem -t 4 -L 80 -M -C ref.fasta adaptor_MergeTrim.fa > s_align.sam

samtools view s_align.sam |  awk 'BEGIN{ OFS= "\t" }{ for (i=12; i<=NF; i++) { if ($i ~ /^XJ/) { split($i,a,":"); print a[3],$3 }} }'| grep -v '*'| sort  > s_align.txt;

zcat barcord_insert_obs.tsv.gz |  python2.7 checkMajority_pipe.py -m 3 | awk 'BEGIN{ OFS="\t" }{ print $3,$2,$1}' | gzip -c > assignment.m3f05.tsv.gz 


2. Count 

paste <(zcat ${ID}*R1_001.fastq.gz) <(zcat .${ID}*R3_001.fastq.gz) <(zcat ${ID}*R2_001.fastq.gz) | awk '{ count+=1; if ((count == 1) || (count == 3)) { print $1 } else { print $1"$$BC$$"$2$3 }; if (count == 4) { count=0 } }' > ${ID}.fastq &

python2.7 SplitFastQdoubleIndexBAM.py -s NUM1 -l NUM2 -m NUM3 -i ../index.lst --remove --summary ../${SAMPLE}.fastq
python2.7MergeTrimReadsBAM.py --mergeoverlap -f FORWARD-ADAPTOR -s REVERSE-ADAPTOR demultiplex.bam 

samtools view -F -r MergeTrim.bam | awk 'BEGIN{ OFS= "\t" }{ for (i=12; i<=NF; i++) { if ($i ~ /^XJ:Z:/) print $10,substr($i,6) }}' | sort | gzip -c  > raw_counts/${SAMPLE}.tsv.gz
zcat raw_counts/${SAMPLE}.tsv.gz  | grep -v "N" | awk 'BEGIN{ OFS="\t" }{ if (length($1) == 15) { print } }' | gzip -c > filter_counts/${SAMPLE}.filtered.tsv.gz

#############
for i in filtered_counts/*.tsv.gz; do echo $(basename $i); join -1 1 -2 1 -t"$(echo -e '\t')" <( zcat $i | cut -f 1 | sort | uniq -c | awk 'BEGIN{ OFS="\t" }{ print $2,$1 }' ) <( zcat assignment/assignment.m3f05.tsv.gz | awk 'BEGIN{ OFS="\t" }{ print $2,$1,$3 }' | sort ) | gzip -c > result/$(basename $i | cut -f 1 -d'.' | awk '{ split($1,a,"-"); print a[1]"_"a[3]"-"a[2] }').tsv.gz; done