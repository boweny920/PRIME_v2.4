#!/usr/bin/env python

import pandas as pd 
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--lims_info_table', help="lims_info.csv fetched from the LIMS api")
parser.add_argument('--RoboIndex_samplesheet')
args=parser.parse_args()

df_lims_info = pd.read_csv(args.lims_info_table)
df_samplesheet = df_lims_info[["libID", "genomeVersion", "speciesName", "orderType","resultPaths", "requestingDepartment", "prnOrderNo", "requester"]].fillna("NA")
df_samplesheet_merged = df_samplesheet.groupby(df_samplesheet.columns.to_list(), as_index=False).first()
df_samplesheet_merged["requestingDepartment"] = df_samplesheet_merged["requestingDepartment"].str.replace("Lab","").str.rstrip().str.replace(" ","_")
# df_samplesheet_merged["speciesName"] = df_samplesheet_merged["speciesName"].str.replace(" ","_") # This is to make speciesName on lims have spaces within
# print([v.replace(" ","_").split("_(")[0] for v in df_samplesheet_merged["speciesName"].to_list()])
df_samplesheet_merged["speciesName"] = [v.replace(" ","_").split("_(")[0] for v in df_samplesheet_merged["speciesName"].to_list()]

#Adding genome annotation version from RoboIndex SampleSheet
available_genomesAnnotations={}
df_robosamplesheet = pd.read_csv(args.RoboIndex_samplesheet)
for index, row in df_robosamplesheet.iterrows():
     available_genomesAnnotations.setdefault( str(row["id"]), str(row["annotation_version"]) )
annotation_list =[]
for genomeversion in df_samplesheet_merged["genomeVersion"].to_list():
    if genomeversion in available_genomesAnnotations:
        annotation_list.append(available_genomesAnnotations[genomeversion])
    else:
        # sys.exit(f"Genome reference {genomeversion} not in Reference Index collection.")
        annotation_list.append("NA")
        
df_samplesheet_merged["annotation"] = annotation_list

df_samplesheet_merged.to_csv("samplesheet.csv", index=False)
