#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -M jdg323@drexel.edu
#$ -P nsftuesPrj
#$ -q all.q@@amdhosts

# ---- Keep the following
. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
# ---- Keep the foregoing

ids=( \
# SRR492065 \
SRR492066 \
SRR492182 \
)
for i in "${ids[@]}"
do :
	echo $i
	fastq-dump -O data -A $i
done
