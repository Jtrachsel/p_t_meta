#!/bin/bash

#SBATCH --job-name="SAMPLE"                            # name of the job submitted
#SBATCH -p brief-low                                   # name of the queue you are submitting to
#SBATCH -n 30                                          # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -N 1                                           # number of nodes
#SBATCH --mem=300G                                     # memory allocation 
#SBATCH -t 2:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH --mail-user=YOUREMAILHERE@email.com            # will receive an email when job starts, ends or fails
#SBATCH --mail-type=BEGIN,END,FAIL                     # will receive an email when job starts, ends or fails
#SBATCH -o "stdout.%j.%N"                              # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                              # optional but it prints our standard error

# ENTER COMMANDS HERE:

module load java/1.8.0_121

set -e # I think this means if anything fails it immediately exits and doesnt keep plowing ahead.

#Written by Brian Bushnell
#Last updated March 4, 2019

#This script is designed to preprocess data for assembly of overlapping 2x150bp reads from Illumina HiSeq 2500.
#Some numbers and steps may need adjustment for different data types or file paths.
#For large genomes, tadpole and bbmerge (during the "Merge" phase) may need the flag "prefilter=2" to avoid running out of memory.
#"prefilter" makes these take twice as long though so don't use it if you have enough memory.
#The "rm SAMPLE_temp.fq.gz; ln -s reads.fq.gz SAMPLE_temp.fq.gz" is not necessary but added so that any pipeline stage can be easily disabled,
#without affecting the input file name of the next stage.

#interleave reads
reformat.sh in=SAMPLE_R#.fastq.gz out=SAMPLE.fq.gz

#Link the interleaved input file as "SAMPLE_temp.fq.gz"
ln -s SAMPLE.fq.gz SAMPLE_temp.fq.gz

# --- Preprocessing ---

#Remove optical duplicates
clumpify.sh in=SAMPLE_temp.fq.gz out=SAMPLE.clumped.fq.gz dedupe optical
rm SAMPLE_temp.fq.gz; ln -s SAMPLE.clumped.fq.gz SAMPLE_temp.fq.gz

#Remove low-quality regions
filterbytile.sh in=SAMPLE_temp.fq.gz out=SAMPLE.filtered_by_tile.fq.gz
rm SAMPLE_temp.fq.gz; ln -s SAMPLE.filtered_by_tile.fq.gz SAMPLE_temp.fq.gz

#Trim adapters.  Optionally, reads with Ns can be discarded by adding "maxns=0" and reads with really low average quality can be discarded with "maq=8".
bbduk.sh in=SAMPLE_temp.fq.gz out=SAMPLE.trimmed.fq.gz ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=adapters ftm=5 ordered
rm SAMPLE_temp.fq.gz; ln -s SAMPLE.trimmed.fq.gz SAMPLE_temp.fq.gz

#Remove synthetic artifacts and spike-ins by kmer-matching.
bbduk.sh in=SAMPLE_temp.fq.gz out=SAMPLE.filtered.fq.gz k=31 ref=artifacts,phix ordered cardinality
rm SAMPLE_temp.fq.gz; ln -s SAMPLE.filtered.fq.gz SAMPLE_temp.fq.gz
