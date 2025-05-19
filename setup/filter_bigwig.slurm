#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=16gb
#SBATCH --time=8:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH --array=1-6
#SBATCH --job-name=bigWig
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o ~/log/bigwig.%A.%a.out
#SBATCH -e ~/log/bigwig.%A.%a.err

source /programs/biogrids.shrc
source $HOME/miniconda3/etc/profile.d/conda.sh

conda activate ucsc-tools

cd $SLURM_SUBMIT_DIR

STIMLIST=./data/stims.txt
STIM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID' $STIMLIST )
INPUT=./data/bigwig/${STIM}.mRp.clN.bigWig
CHROMSIZES=./data/chrom.sizes
TEMP=./data/${STIM}.bedgraph
TEMPFILT=./data/${STIM}_filtered.bedgraph
OUTPUT=../inst/extdata/${STIM}_filtered.bigWig

bigWigToBedGraph $INPUT $TEMP

awk '$4 > 0.1' $TEMP > $TEMPFILT

bedGraphToBigWig $TEMPFILT $CHROMSIZES $OUTPUT

rm $TEMP $TEMPFILT

conda deactivate
