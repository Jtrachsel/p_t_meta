#!/bin/bash

#SBATCH --job-name="ALL_QC"                            # name of the job submitted
#SBATCH -p mem                                         # name of the queue you are submitting to
#SBATCH -n 40                                          # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -N 1                                           # nodes
#SBATCH --mem=1000G                                    # memory allocation  
#SBATCH -t 72:00:00                                    # time allocated for this job hours:mins:seconds
#SBATCH --mail-user=YOUR_EMAIL@email.com               # will receive an email when job starts, ends or fails
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
#The "rm ALL_temp.fq.gz; ln -s reads.fq.gz ALL_temp.fq.gz" is not necessary but added so that any pipeline stage can be easily disabled,
#without affecting the input file name of the next stage.

#interleave reads
#reformat.sh in=ALL_R#.fastq.gz out=ALL.fq.gz

# 92 fecal shotgun metagenomic libraries
# sequenced on HiSeq3000 2x150 PE
# I think avg insert sizes is ~280?

# this takes the R1 and R2 files and creates interleaved fastq.gz files for each
for x in *_R1.fastq.gz 
do
#echo "${x%_R1*}_R2.fastq.gz"
reformat.sh in1="$x" in2="${x%_R1*}_R2.fastq.gz" out="${x%_R1*}.fq.gz"
done

# combines all interleaved reads into one file
cat *.fq.gz > ALL.fq.gz

#Link the interleaved input file as "ALL_temp.fq.gz"
ln -s ALL.fq.gz ALL_temp.fq.gz

# --- Preprocessing ---

#Remove optical duplicates
clumpify.sh in=ALL_temp.fq.gz out=ALL.clumped.fq.gz dedupe optical
rm ALL_temp.fq.gz; ln -s ALL.clumped.fq.gz ALL_temp.fq.gz

#Remove low-quality regions
filterbytile.sh in=ALL_temp.fq.gz out=ALL.filtered_by_tile.fq.gz
rm ALL_temp.fq.gz; ln -s ALL.filtered_by_tile.fq.gz ALL_temp.fq.gz

#Trim adapters.  Optionally, reads with Ns can be discarded by adding "maxns=0" and reads with really low average quality can be discarded with "maq=8".
bbduk.sh in=ALL_temp.fq.gz out=ALL.trimmed.fq.gz ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=adapters ftm=5 ordered
rm ALL_temp.fq.gz; ln -s ALL.trimmed.fq.gz ALL_temp.fq.gz

#Remove synthetic artifacts and spike-ins by kmer-matching.
bbduk.sh in=ALL_temp.fq.gz out=ALL.filtered.fq.gz k=31 ref=artifacts,phix ordered cardinality
rm ALL_temp.fq.gz; ln -s ALL.filtered.fq.gz ALL_temp.fq.gz

#Decontamination by mapping can be done here.
#JGI removes these in two phases:
#1) common microbial contaminants (E.coli, Pseudomonas, Delftia, others)
#2) common animal contaminants (Human, cat, dog, mouse)

bbsplit.sh in=SAMPLE_temp.fq.gz ref=turkey.fa basename=SAMPLE_out_%.fq.gz outu=SAMPLE_clean.fq.gz int=t
rm SAMPLE_temp.fq.gz; ln -s SAMPLE_clean.fq.gz SAMPLE_temp.fq.gz


#Error-correct phase 1
bbmerge.sh in=ALL_temp.fq.gz out=ALL.ecco.fq.gz ecco mix vstrict ordered ihist=ihist_merge1.txt
rm ALL_temp.fq.gz; ln -s ALL.ecco.fq.gz ALL_temp.fq.gz

#Error-correct phase 2
clumpify.sh in=ALL_temp.fq.gz out=ALL.eccc.fq.gz ecc passes=4 reorder
rm ALL_temp.fq.gz; ln -s ALL.eccc.fq.gz ALL_temp.fq.gz

#Error-correct phase 3
#Low-depth reads can be discarded here with the "tossjunk", "tossdepth", or "tossuncorrectable" flags.
#For very large datasets, "prefilter=1" or "prefilter=2" can be added to conserve memory.
#Alternatively, bbcms.sh can be used if Tadpole still runs out of memory.
tadpole.sh in=ALL_temp.fq.gz out=ALL.ecct.fq.gz ecc k=62 ordered tossjunk prefilter=1
rm ALL_temp.fq.gz; ln -s ALL.ecct.fq.gz ALL_temp.fq.gz

#Normalize
#This phase can be very beneficial for data with uneven coverage like metagenomes, MDA-amplified single cells, and RNA-seq, but is not usually recommended for isolate DNA.
#So normally, this stage should be commented out, as it is here.
bbnorm.sh in=ALL_temp.fq.gz out=normalized.fq.gz target=100 hist=khist.txt peaks=peaks.txt
rm ALL_temp.fq.gz; ln -s normalized.fq.gz ALL_temp.fq.gz


