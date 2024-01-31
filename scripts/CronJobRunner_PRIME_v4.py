
import logging
import os
import requests
from subprocess import check_output, CalledProcessError, STDOUT, call
import argparse
import pandas as pd

parser=argparse.ArgumentParser(description="CrobJob to start the PRIME pipeline")
parser.add_argument('--samplesheet',
        help='Path to samplesheet used to ROBOINDEX update')
args = parser.parse_args()

def system_call(command: list) -> tuple:
    try:
        output = check_output(command, stderr=STDOUT).decode()
        success = True
    except CalledProcessError as e:
        output = e.output.decode()
        success = False
    return output, success

NGS_LIMS = 'https://lims.stowers.org/zanmodules/molecular-biology/ngs'
API_TOKEN = 'ca7952666a03dd4e59d0cd59e39fecc7'
requests.packages.urllib3.util.ssl_.DEFAULT_CIPHERS = 'ALL:@SECLEVEL=1'

def fetch_lims(fc_id: str,
                ngs_lims: str = NGS_LIMS,
                api_token: str = API_TOKEN) -> dict:
    """
    Get the run info from LIMS, mostly just need the MOLNG ID for Miseq runs...
    """
    # Get run info from lims
    ngs_lims = NGS_LIMS
    api_token = API_TOKEN  # User specific
    header = {'x-zan-apitoken': f'{api_token}', 'Accept': 'application/json'}
    run_info = requests.get(f'{ngs_lims}/flowcells/{fc_id}/samples',
                            headers=header, verify=False)
    if not run_info.ok:
        print(f"Issue fetching {fc_id} in LIMS")
        return None
    lims_data = run_info.json()  # Close request
    unique_molng = list(set([x['prnOrderNo'] for x in lims_data['samples']]))
    unique_id = list(set([x['orderID'] for x in lims_data['samples']]))
    unique_orderType = list(set([x['orderType'] for x in lims_data['samples']]))
    unique_analysisPath = list(set([x['resultPaths'] for x in lims_data['samples']]))
    unique_genome_ver = list(set([x['genomeVersion'] for x in lims_data['samples']]))
    # unique_species = list(set([x['speciesName'].replace(" ","_").split("_(")[0] for x in lims_data['samples']]))
       
    lims_info = {
        'FCID': lims_data['FCID'],
        'limsID': lims_data['prnId'],
        'MOLNGs': unique_molng[0],
        'orderID': unique_id[0],
        'orderType' : unique_orderType[0],
        'nanalysis_path' : unique_analysisPath[0],
        'genome_ver' : unique_genome_ver[0],
        # 'species' : unique_species[0] 
    }
    return lims_info

def get_fc_id(fc_path: str) -> str:
    """
    Get flowcell ID from flowcell path
    """
    fc_id = fc_path.split('_')[-1]
    return fc_id

def check_log(fc_dict: dict, logfile):
    """
    Check log file if fc id is already there, don't add if it is
    """
    seen_key = []
    with open(logfile, 'r') as f:
        for line in f:
            for fcid in fc_dict.keys():
                if fcid in line:
                    seen_key.append(fcid)
    fc_dict = {key: val for key, val in fc_dict.items() if key not in seen_key}
    return fc_dict

def new_fc_pipeline_run(logfile, copy_complete_txt_name, data_folder, machine_type, available_genomes):

    logging.basicConfig(filename=logfile,
                        level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S')

    output, success = system_call([
            "/n/ngs/tools/primary_smk/bin/src/fd-v8.2.1-i686-unknown-linux-musl/fd", # copy the fd to one of PRIME's folders
            copy_complete_txt_name, data_folder, 
            "-d", "2",
            "--changed-within", "128h", "--threads", "1", 
            "--type", "f"
        ])
    if not success:
            print("Something isn't working! Check the code")

    paths = [os.path.dirname(out) for out in output.split("\n") if os.path.dirname(out) != '']
    fc_dict = {}
    for path in paths:
        fc_id = get_fc_id(path)
        lims_info = fetch_lims(fc_id)
        if lims_info != None:
            if "Targeted Deep" not in lims_info["orderType"]: # Making sure it is not genome engineering orders!! 
                # fc_dict[fc_id] = {"fc_path":path, "orderType":lims_info["orderType"], "nanalysis_path":lims_info['nanalysis_path'], "reference_genome":lims_info["species"]+"_"+lims_info["genome_ver"]}
                fc_dict[fc_id] = {"fc_path":path, "orderType":lims_info["orderType"], "nanalysis_path":lims_info['nanalysis_path'], "reference_genome": lims_info["genome_ver"]} # Only using the genome version as a checking key for now
    
    fc_dict = check_log(fc_dict, logfile) # Remove fc from dict if it's alraedy in the log file
    
    if len(fc_dict) > 0:
        for fc in fc_dict:
            logging.info(f"Discovered {fc} - fc_path: {fc_dict[fc]['fc_path']}")
            os.chdir(fc_dict[fc]['fc_path'])
            
            if fc_dict[fc]["reference_genome"] in available_genomes:
                only_fastq = ""
            else: 
                only_fastq = "--only_cp_fastq True"
            
            if "10x_scRNA-Seq_NextGEM" in fc_dict[fc]["orderType"]:
            # if it is 10X single cell:
                cmd = f"source ~/.bash_profile; source /etc/profile; source ~/.bashrc; PATH=\"/usr/local/bin:$PATH\"; nextflow run /n/ngs/tools/PRIME/RUN/main.nf \
                    --fcid {fc} --fcpath {fc_dict[fc]['fc_path']} --machine_type {machine_type} --run_type 10X -profile slurm {only_fastq} > \
                    /n/ngs/tools/PRIME/logs/machine_run_logs/{fc}.{machine_type}.nf.log"
            else:
                cmd = f"source ~/.bash_profile; source /etc/profile; source ~/.bashrc; PATH=\"/usr/local/bin:$PATH\"; nextflow run /n/ngs/tools/PRIME/RUN/main.nf \
                    --fcid {fc} --fcpath {fc_dict[fc]['fc_path']} --machine_type {machine_type} --run_type RNA-Seq -profile slurm {only_fastq} > \
                    /n/ngs/tools/PRIME/logs/machine_run_logs/{fc}.{machine_type}.nf.log"
            
            start_email = f"mail -s 'Flowcell: {fc}; Order Type: {fc_dict[fc]['orderType']}; Machine: {machine_type}; PRIME execution STARTED' -r by2747@stowers.org by2747@stowers.org,mpe@stowers.org,hhassan@stowers.org,mcm@stowers.org" #,mpe@stowers.org
            call(start_email, shell=True)
            print(cmd)
            return_code = call(cmd, shell=True)
            if return_code == 0:
                logging.info(f"PRIME Nextflow Pipeline executed successfully {fc}")
                print("Command executed successfully")
                # email = f"echo {fc} fastqs stored at : https://webfs/{fc_dict[fc][1]} | mail -s '{fc} Genome Engineering Pipeline execution COMPLETED' -r by2747@stowers.org by2747@stowers.org,HHassan@stowers.org,kjw@stowers.org"
            else:
                print("Command failed with exit code", return_code)
               
            
    else:
        logging.info("No new runs")

PRIME_log = "/n/ngs/tools/PRIME/logs/PRIME_Orders.log"

## Make dictionary of availabe reference genomes, thier versions, and annotation versions in current SIMR collection
available_genomes = {}
df_samplesheet = pd.read_csv(args.samplesheet)
for index, row in df_samplesheet.iterrows():
    # available_genomes.setdefault( str(row["name"]).replace(" ","_") + "_" + str(row["id"]), row["annotation_version"] )
    available_genomes.setdefault( str(row["id"]), row["annotation_version"] ) # Use only genmome version as key for now
    
## Check NextSeq2000 runs ###
machine_type = "NextSeq2K"
copy_complete_txt_name = "CopyComplete.txt"
data_folder = "/n/ngs/data/NextSeq2K"
new_fc_pipeline_run(PRIME_log, copy_complete_txt_name, data_folder, machine_type, available_genomes) 

### Check G4 runs ###
machine_type = "G4"
copy_complete_txt_name = "transfer_complete.txt"
data_folder = "/n/ngs/data/G4"
new_fc_pipeline_run(PRIME_log, copy_complete_txt_name, data_folder, machine_type, available_genomes) 

### Check Aviti runs ###
machine_type = "Aviti"
copy_complete_txt_name = "RunUploaded.json" 
data_folder = "/n/ngs/data/Aviti/AV230402"
new_fc_pipeline_run(PRIME_log, copy_complete_txt_name, data_folder, machine_type, available_genomes) 