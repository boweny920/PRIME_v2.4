#!/usr/bin/env python
import sys
import pandas as pd 
import os 
import argparse
import requests

parser = argparse.ArgumentParser()
parser.add_argument('--output_dir', default="./", help="Output directory, default=cwd, default='./'")
parser.add_argument('--lims_info_csv', help='provide fcid')
parser.add_argument('--libID', help='library ID')
parser.add_argument('--fastqDir', help='Path to fastqs')
args=parser.parse_args()

df_lims_info = pd.read_csv(args.lims_info_csv)
df_lims_info = df_lims_info[df_lims_info["libID"] == args.libID][[ "libID","analysisGoals","genomeVersion"]].drop_duplicates()

## Get the sample info ##
analysis_goal_string = df_lims_info[df_lims_info["libID"] == args.libID]["analysisGoals"].to_list()[0]
config_csv_ready_sampleInfo = "\n".join([v.split("Pool_")[0].replace(" ","") for v in analysis_goal_string.split("\n") if "BC00" in v and str(args.libID) in v])
# config_csv_ready_sampleInfo e.g.
# S62457,BC001,L63060,
# S62458,BC002,L63060,
# S62459,BC003,L63060,
# S62460,BC004,L63060,

## Getting the 10x reference info and probe file, current probe files only support Ens_98 mm10 & hg38 ##
if len(df_lims_info[df_lims_info["libID"] == args.libID]["genomeVersion"].to_list()) > 1:
    sys.exit(f"More than 1 reference for {args.libID}!")
if "hg" in df_lims_info[df_lims_info["libID"] == args.libID]["genomeVersion"].to_list()[0]:
    reference = "/n/analysis/indexes/hg38/annotation/Ens_98/10x/cellranger-6.0.1/"
    probe_file = "/n/core/Bioinformatics/analysis/CompBio/boweny/nf-Pipeline/nextflow-GenomeIndex_Ver_4/bin/external/tenx_feature_references/targeted_panels/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"
elif "mm" in df_lims_info[df_lims_info["libID"] == args.libID]["genomeVersion"].to_list()[0]:
    reference = "/n/analysis/indexes/mm10/annotation/Ens_98/10x/cellranger-6.0.1"
    probe_file = "/n/core/Bioinformatics/analysis/CompBio/boweny/nf-Pipeline/nextflow-GenomeIndex_Ver_4/bin/external/tenx_feature_references/targeted_panels/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv"
else:
    sys.exit("Currently only mm10 and hg38 are supported with 10 FLEX, none detected in lims_info.csv")

config_csv = f"[gene-expression]\n\
reference,{reference},,,,,\n\
probe-set,{probe_file},,,,,\n\
,,,,,,\n\
[libraries]\n\
fastq_id,fastqs,lanes,feature_types,subsample_rate\n\
{args.libID},{args.fastqDir},,Gene Expression,\n\
,,,,,,\n\
[samples]\n\
sample_id,probe_barcode_ids,description,expect_cells,force_cells\n\
{config_csv_ready_sampleInfo}\n\
,,,,,,"

with open(f"{args.libID}.config.csv","w") as newfile:
    newfile.write(config_csv)
    