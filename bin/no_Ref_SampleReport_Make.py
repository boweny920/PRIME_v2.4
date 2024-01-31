#!/usr/bin/env python 

import pandas as pd 
import argparse
import sys
import os

parser = argparse.ArgumentParser()
parser.add_argument('--output_dir', default="./", help="Output directory, default=cwd, default='./'")
parser.add_argument('--driver_csv', help='provide fcid')
args=parser.parse_args()

df_driver = pd.read_csv(args.driver_csv)

df_sampleReport = df_driver[["Order","OrderType","SampleName","LibraryID","IndexSequences","Read","Reference","Species","Lab","Type","ReadLength"]]
df_sampleReport["Output"] = df_driver["fileNames"]
df_sampleReport["FastqPath"] = df_driver["resultPaths"].astype(str) + "/" + df_driver["fileNames"].astype(str) + ".fastq.gz"
df_sampleReport["TotalReads"] = ["NA-NotCalculated"]*len(df_sampleReport)
df_sampleReport["AlignPercent"] = ["NA-NotCalculated"]*len(df_sampleReport)
# df_sampleReport[["IndexSequence1", "IndexSequence2"]] = df_sampleReport["IndexSequences"].str.split(pat="-",expand=True)
df_sampleReport["IndexSequence1"], df_sampleReport["IndexSequence2"] = zip(*df_sampleReport["IndexSequences"].str.split("-").apply(lambda x: x + [None] if len(x) == 1 else x).tolist())
# Handling cases where no "-" is present
df_sampleReport["IndexSequence2"].fillna("NA", inplace=True)

col_order = ["Output","Order","OrderType","SampleName","LibraryID","IndexSequence1","IndexSequence2","Read","Reference","Lab","TotalReads","AlignPercent","Type","ReadLength","Species","FastqPath"]
df_sampleReport = df_sampleReport[col_order]
df_sampleReport.sort_values(by='Output', inplace=True)
df_sampleReport.to_csv("Prime_Sample_Report.csv", index=False)

for molng_id in df_sampleReport["Order"].unique():
    df_molng = df_sampleReport[df_sampleReport["Order"]==molng_id]
    nanalysis_path =df_driver[df_driver['Order']==molng_id]['resultPaths'].unique()
    
    if len(nanalysis_path)>1:
        sys.exit(f"more than one nanalysis_path associated to one {molng_id}")
    nanalysis_path = nanalysis_path[0]
    if os.path.exists(nanalysis_path) == False:
        os.makedirs(nanalysis_path)
    
    df_molng.to_csv(os.path.join(nanalysis_path,"Sample_Report.csv"), index=False)
