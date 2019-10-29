#!/bin/bash

#SBATCH --job-name="SAMPLE"                             # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -n 8                                            # number of cores/tasks 
#SBATCH -N 1
#SBATCH --mem=100G                                      # memory allocation
#SBATCH -t 47:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH --mail-user=YOUREMAIL@email.com                 # will receive an email when job fails
#SBATCH --mail-type=FAIL                                # will receive an email when job fails
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error

# ENTER COMMANDS HERE:

# activate the conda environment
source ~/iREP/bin/activate

gunzip SAMPLE_mapped.sam.gz

# allows iRep to tolerate 7 mismatches, over ~150bp reads this is about 95%
iRep -f *.fa -s SAMPLE_mapped.sam -o SAMPLE.out -mm 7 -t 8 --no-plot --sort

