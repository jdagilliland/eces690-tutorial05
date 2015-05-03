#!/bin/bash
# tell SGE to use bash for this script
#$ -S /bin/bash
# execute the job from the current working directory, i.e. the directory in which the qsub command is given
#$ -cwd
# set email address for sending job status
#$ -M jdg323@drexel.edu
# project - basically, your research group name with "Grp" replaced by "Prj"
#$ -P nsftuesPrj
# select parallel environment, and number of job slots
#$ -pe openmpi_ib 16
# request 15 min of wall clock time "h_rt" = "hard real time" (format is HH:MM:SS, or integer seconds)
#$ -l h_rt=24:00:00
# a hard limit 8 GB of memory per slot - if the job grows beyond this, the job is killed
#$ -l h_vmem=8G
##$ -pe shm 32-64 #for parallel processing
# want at least 6 GB of free memory
#$ -l mem_free=6G
# select the queue all.q, using hostgroup @intelhosts
#$ -q all.q@@amdhosts

PATH=/mnt/HA/groups/nsftuesGrp/.local/bin:$PATH

# ---- Keep the following
. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
# ---- Keep the foregoing

# Set up variables to keep organized.
# DATADIR=/mnt/HA/groups/nsftuesGrp/data/tutorial5
DATADIR=$HOME/devo/tutorial05/data
INDIR=${DATADIR}/fastq-reads
mkdir -p $INDIR
OUTDIR=${TMP}/idbaout
mkdir -p $OUTDIR

# Download the requisite data.
ids=( \
SRR492065 \
SRR492066 \
SRR492182 \
)
# for i in "${ids[@]}"
# do :
# 	echo $i
# 	fastq-dump -O $INDIR -A $i
# done

# Assemble the reads *de novo*.
idba_ud \
	-l $INDIR/SRR492065.fastq \
	$INDIR/SRR492066.fastq \
	$INDIR/SRR492185.fastq \
	-o $OUTDIR

# Clean up, clean up, everybody clean up.
# rm -rf $INDIR
mv $OUTDIR/* $DATADIR/
