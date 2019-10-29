#!/bin/bash

#SBATCH --job-name="SAMPLE"                            # name of the job submitted
#SBATCH -p mem-low                                     # name of the queue you are submitting to
#SBATCH -n 40                                          # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -N 1                                           # number of nodes
#SBATCH --mem=200G                                     # memory allocation
#SBATCH -t 2:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH --mail-user=youremail@email.com                # will receive an email when job starts, ends or fails
#SBATCH --mail-type=BEGIN,END,FAIL                     # will receive an email when job starts, ends or fails
#SBATCH -o "stdout.%j.%N"                              # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                              # optional but it prints our standard error

# ENTER COMMANDS HERE:

module load java/1.8.0_121
module load samtools
module load pigz

set -e # I think this means if anything fails it immediately exits and doesnt keep plowing ahead.



# FIRST GENERATING AN INDEXED REFERENCE WITH 
# bbmap.sh ref=metaspades_500_scaff.fa

# then all jobs will just use this pre-indexed reference. I am doing this with an interactive shell prior to submitting the SLURMS

# not specifying reference because of preindexing listed above.


bbmap.sh in=../hiseq/SAMPLE.filtered.fq.gz outm=SAMPLE_mapped.sam.gz bs=SAMPLE.bamscript qtrim=10 untrim ambig=all pigz unpigz

# this bamscript will generate sorted & indexed bam files
chmod u+x SAMPLE.bamscript
./SAMPLE.bamscript
