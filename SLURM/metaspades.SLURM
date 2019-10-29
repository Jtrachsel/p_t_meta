#!/bin/bash

#SBATCH --job-name="metaspades"                         # name of the job submitted
#SBATCH -p mem                                          # name of the queue you are submitting to
#SBATCH -n 100                                          # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -N 1                                            # number of nodes
#SBATCH --mem=1400G                                     # memory allocation
#SBATCH -t 167:00:00                                    # time allocated for this job hours:mins:seconds
#SBATCH --mail-user=YOUREMAILHERE@email.com             # will receive an email when job starts, ends or fails
#SBATCH --mail-type=BEGIN,END,FAIL                      # will receive an email when job starts, ends or fails
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error


module load java/1.8.0_121

# convert from interleaved to R1 R2 files
reformat.sh in=ALL_normalized.fq.gz out1=ALL_normalized_1.fq.gz out2=ALL_normalized_2.fq.gz

metaspades.py -o run1 -1 ALL_normalized_1.fq.gz -2 ALL_normalized_2.fq.gz --pacbio ALL_pacbio.fq.gz --only-assembler -m 1400 --tmp-dir "$TMPDIR" -t 100 -k 25,55,95,125
