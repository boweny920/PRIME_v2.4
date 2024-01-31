import argparse
import os 
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fastq')
args=parser.parse_args()

new_name_pre = str(args.fastq).replace(".fastq.gz","").split("_")
libID=new_name_pre[0]
read=new_name_pre[1]

newfastq = libID+"_S1_L001_"+ read + "_001.fastq.gz"
#print(newfastq)

shutil.move(args.fastq, newfastq)
