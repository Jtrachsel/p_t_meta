#!/bin/bash

#SBATCH --job-name="7_compile"                          # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -n 8                                            # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -N 1                                            # nodes
#SBATCH --mem=32G                                       # memory allocation
#SBATCH -t 47:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH --mail-user=yourmail@mail.com                   # will receive an email when job fails
#SBATCH --mail-type=FAIL                                # will receive an email when job fails
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error

# ENTER COMMANDS HERE:

source ~/iREP/bin/activate


# iRep -f BIN.fa -s BIN*sam -o BIN.out -mm 7 -t 8 --noplot




iRep_filter.py -c 3 -w 93 -t *tsv --long > ALL_RESULTS_long
