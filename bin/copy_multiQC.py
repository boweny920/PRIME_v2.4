#!/usr/bin/env python 

import pandas as pd 
import argparse
import os 
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--driver_csv')
parser.add_argument('-m', '--multiqc')
args=parser.parse_args()
#use for test: /n/ngs/data/NextSeq2K/230919_VH00629_133_AAC7NV3HV/Sample_Report.csv

df = pd.read_csv(args.driver_csv)
for nanalysis_path in df["resultPaths"].unique():
    if os.path.exists(nanalysis_path) == False:
        os.makedirs(nanalysis_path)
    shutil.copy2(args.multiqc, nanalysis_path)
