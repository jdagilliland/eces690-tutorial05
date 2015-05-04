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

# Load bowtie2 module.
module load bowtie2/2.2.5
# ---- Keep the foregoing

# Set up variables to keep organized.
# DATADIR="$HOME/devo/tutorial05/data/idbaout"
DATADIR="$HOME/devo/tutorial05/data"
INDEXNAME="gut-biome"
INDEXDIR="$TMP/index"
mkdir -p $INDEXDIR

# Build the index.
bowtie2-build "$DATADIR/contig-100.fa","$DATADIR/contig-20.fa","$DATADIR/contig-40.fa","$DATADIR/contig-60.fa","$DATADIR/contig-80.fa" "$INDEXDIR/$INDEXNAME"

ids=( \
SRR492065 \
SRR492066 \
SRR492182 \
)

for i in "${ids[@]}"
do :
	echo $i
	bowtie2 -x "$DATADIR/$INDEXNAME" -U "${DATADIR}/${i}.fastq" \
	-S "${INDEXDIR}/${i}.sam"
	# bowtie2 -x "$INDEXDIR/$INDEXNAME" -U "${DATADIR}/${i}.fastq" \
	# -S "${INDEXDIR}/${i}.sam"
done

# Clean up.
mv $INDEXDIR/* $DATADIR/
rmdir $INDEXDIR
