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
PATH=/mnt/HA/groups/nsftuesGrp/.local/bin:$PATH

# Load bowtie2 module.
# module load bowtie2/2.2.5
# ---- Keep the foregoing

echo $PATH
idba_ud --help

# ids=( \
# SRR492065 \
# # SRR492066 \
# # SRR492182 \
# )
