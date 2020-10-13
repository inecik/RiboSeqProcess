"""
Analysis of Paired End Sequencing of Disome Profiling
"""


import os
import subprocess
import sys
import shutil
import pandas as pd
from datetime import datetime
import multiprocessing

# Welcoming Message
print("Analysis of Paired End Sequencing of Disome Profiling")

# Excel file for annotated supplied by the user
excel_path = "examples/example_annotation.xlsx"  # todo: sys.argv[0]
raw_annotated = pd.read_excel(excel_path)  # Read the excel file as dataframe object


# Create and set working directory
parent_directory = os.path.dirname(excel_path)
output_dir = os.path.join(parent_directory, "output")  # output directory is next to the excel file provided
output_dir += datetime.now().strftime("-%Y.%m.%d-%H.%M.%S")  # Timestamp is added for data security and convenience
os.mkdir(output_dir)  # Create such directory
os.chdir(output_dir)  # Working environment is set to the newly created directory

a=subprocess.run("which python", shell=True)

# Create necessary folders in working directory
create_directories = ["01_Combined_Data",
                      "02_Cutadapt_Trimmed_Data",
                      "03_Bowtie2_rRNA_Depleted_Data",
                      "04_STAR_Alignment_Data",
                      "05_FastQC_Analysis",
                      "06_Reads_Assignment"]
create_directories.extend([os.path.join("04_STAR_Alignment_Data", i) for i in pd.unique(raw_annotated["Repeat"])])
create_directories.extend([os.path.join("06_Reads_Assignment", i) for i in pd.unique(raw_annotated["Repeat"])])
for i in create_directories:
    os.mkdir(os.path.join(os.getcwd(), i))


# Trim everything else other than footprint. The structure of the final library is as follows.
# [Read1 sequencing primer annealing seq] - NN - footprint - NNNNN - BBBBB - [Read2 sequencing primer annealing seq]
for excel_rep in pd.unique(raw_annotated["Repeat"]).tolist():
    excel_read_s = raw_annotated[raw_annotated["Repeat"] == excel_rep]["Read"].tolist()
    excel_path_s = raw_annotated[raw_annotated["Repeat"] == excel_rep]["Full Path"].tolist()
    trim_output_dir = "02_Cutadapt_Trimmed_Data"
    trim_outputs = [f"{trim_output_dir}/{excel_rep}_{r}_Cutadapt-Trimmed.fastq.gz" for r in excel_read_s]
    # Note always R1 is above R2 in excel file

    # todo: buradan sonraki raporu editleyen scripti anlamadım la


# rRNA removal
for excel_rep in pd.unique(raw_annotated["Repeat"]).tolist():
    excel_read_s = raw_annotated[raw_annotated["Repeat"] == excel_rep]["Read"].tolist()
    excel_path_s = raw_annotated[raw_annotated["Repeat"] == excel_rep]["Full Path"].tolist()
    trim_output_dir = "03_Bowtie2_rRNA_Depleted_Data"
    trim_outputs = [f"{trim_output_dir}/{excel_rep}_{r}_Cutadapt-Trimmed.fastq.gz" for r in excel_read_s]
    cut_adapt_command


    bowtie2 - p8 - t - I
    40 - X
    90 - -no - mixed - x
    '/media/matilde/Matys_Disk/Data_Analysis_new/data_files/indexed_rRNA_23.06.19/indexed_rRNA' - 1
    '02_Cutadapt_Trimmed_Data/'$val
    '_R1_Cutadapt-Trimmed.fastq.gz' - 2
    '02_Cutadapt_Trimmed_Data/'$val
    '_R2_trimmed_umi-aware.fastq.gz' - -un - conc
    '03_Bowtie2_rRNA_Depleted_Data/'$val
    '_norRNA.fastq' - S / dev / null
    2 > '03_Bowtie2_rRNA_Depleted_Data/'$val
    '_Bowtie2_report.txt'



