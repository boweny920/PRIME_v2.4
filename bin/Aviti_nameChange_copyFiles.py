#!/usr/bin/env python

import glob
import argparse
import shutil
import os 

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--fastq_dir')
parser.add_argument('-l', '--libID')
parser.add_argument('-n', '--nanalysis_dir', help='the /n/analysis dir that the fastqs will be copied to')
args=parser.parse_args()

# One library for each sample? According to Mike, yes for now
R1_fastq_file = f"{args.fastq_dir}/{args.libID}_R1.fastq.gz"
R2_fastq_file = f"{args.fastq_dir}/{args.libID}_R2.fastq.gz"

# New file name: MySample_S1_L001_R1_001.fastq.gz [sample name, sample order, lane, read type, and chunk]
newR1_fastq_file = f"{args.fastq_dir}/{args.libID}_S1_L001_R1_001.fastq.gz"
newR2_fastq_file = f"{args.fastq_dir}/{args.libID}_S1_L001_R2_001.fastq.gz"

#Rename to the way cellranger likes it
if os.path.isfile(newR1_fastq_file) == False: 
    shutil.move(R1_fastq_file, newR1_fastq_file)
if os.path.isfile(newR2_fastq_file) == False:
    shutil.move(R2_fastq_file, newR2_fastq_file)

#Copy renamed fastqs to /n/analysis folders
if os.path.exists(args.nanalysis_dir) == False:
    os.makedirs(args.nanalysis_dir)

if os.path.isfile(os.path.join(args.nanalysis_dir, os.path.basename(newR1_fastq_file))) == False:
    shutil.copy2(newR1_fastq_file, args.nanalysis_dir)

if os.path.isfile(os.path.join(args.nanalysis_dir, os.path.basename(newR2_fastq_file))) == False:
    shutil.copy2(newR2_fastq_file, args.nanalysis_dir)
