#!/bin/bash

#SBATCH --job-name="metabat"                            # name of the job submitted
#SBATCH -p mem                                          # name of the queue you are submitting to
#SBATCH -n 40                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -N 1                                            # number of nodes
#SBATCH --mem=125G                                      # memory allocation
#SBATCH -t 47:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH --mail-user=YOUREMAIL@email.com                 # will receive an email when job starts, ends or fails
#SBATCH --mail-type=BEGIN,END,FAIL                      # will receive an email when job starts, ends or fails
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error

# ENTER COMMANDS HERE:

module load metabat
module load checkm/v1.0.11

jgi_summarize_bam_contig_depths --outputDepth depth.txt --pairedContigs paired.txt --minContigLength 1000 --minContigDepth 2 *.bam

metabat -i metaspades_500_scaff.fa -a depth.txt -o run1/bin -v --saveTNF saved.tnf --saveDistance saved.dist
checkm lineage_wf -f run1/CheckM.txt -t 32 -x fa run1/ run1/SCG

# might want to consider doing a different assembly technique... individual assemblies -> bin -> dereplicate
