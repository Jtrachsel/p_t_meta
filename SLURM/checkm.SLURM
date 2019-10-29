#!/bin/bash

#SBATCH --job-name="checkm"                             # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --mem=48G                                       # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 47:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH --mail-user=yourname@mail.com                   # will receive an email when job starts, ends or fails
#SBATCH --mail-type=BEGIN,END,FAIL                      # will receive an email when job starts, ends or fails
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error

# ENTER COMMANDS HERE:

module load checkm

checkm lineage_wf -f out/CheckM.txt -t 16 -x fasta ./ ./out

#
