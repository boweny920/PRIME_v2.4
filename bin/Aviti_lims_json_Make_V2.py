#!/usr/bin/env python

import sys
import pandas as pd 
import os 
import argparse
import requests
import json
import urllib3

requests.packages.urllib3.util.ssl_.DEFAULT_CIPHERS = 'ALL:@SECLEVEL=1'

def fetch_sample_table_from_LIMS(fc_id, molng):
    """
    Get the run info from LIMS. Returns the /n/analysis/ folder where the primary analsysis pipeline copied all the data
    """
    # Get run info from lims
    NGS_LIMS = 'https://lims.stowers.org/zanmodules/molecular-biology/ngs'
    API_TOKEN = 'ca7952666a03dd4e59d0cd59e39fecc7' # Should get a new one for each pipeline, or even each user

    header = {'x-zan-apitoken': f'{API_TOKEN}', 
        'Accept': 'application/json'}
    run_info = requests.get(f'{NGS_LIMS}/flowcells/{fc_id}/samples', headers=header, verify=False)
    if not run_info.ok: 
        sys.exit('Malformed API request. Please double check your flowcell ID')
    lims_data = run_info.json() # Close request
    run_info.close() 
    df_run_info = pd.DataFrame.from_dict(lims_data)
    df_run_info.to_csv("tmp_lims.csv", index=False)
    df_sample_info = pd.concat([df_run_info['readLength'] ,df_run_info['readType'], df_run_info['samples'].apply(pd.Series)], axis=1)
    
    if molng != False:
        df_sample_info = df_sample_info[df_sample_info["prnOrderNo"] == molng]

    return df_sample_info

parser = argparse.ArgumentParser()
# parser.add_argument('-s', '--sample_report', default=None, help= "Input Samplereport.csv associated to order, generated by primary analysis pipeline")
parser.add_argument('-o', '--output_dir', default="./", help="Output directory, default=cwd, default='./'")
parser.add_argument('-l', '--fcid', help='provide fcid')
parser.add_argument("-m", '--molng', default=False, help='provide molng id, this can be just empty if the FCID has only 1 MOLNG-ID')
parser.add_argument("-s", '--samplesheet_Robo', default=False, help='RoboIndex Samplesheet')
args=parser.parse_args()

available_species = {}
df_samplesheet = pd.read_csv(args.samplesheet_Robo)
for index, row in df_samplesheet.iterrows():
    available_species.setdefault( str(row["id"]), str(row["name"]))

df_sample_info = fetch_sample_table_from_LIMS(args.fcid, args.molng)

# Making the correct species name in the downloaded lims_info table
genome_species = []
for v in df_sample_info["genomeVersion"]:
    if v in available_species:
        genome_species.append(available_species[v])
    else:
        genome_species.append("None")
df_sample_info["speciesName"] = genome_species
 
df_sample_info.to_csv("lims_info.csv", index=False)
# requestingDepartment will be pi name 
# requester is the person making the order request