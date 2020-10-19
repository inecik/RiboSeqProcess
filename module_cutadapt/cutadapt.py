"""
Cutadapt for paired end sequencing, runs for both reads
"""


import os
import sys
import subprocess
from multiprocessing import cpu_count
from shutil import which


# CONSTANTS
adapter_seq_read_1 = "ATCGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
adapter_seq_read_2 = None  
# Mati and Kai's protocol does not necessitates to cut read 2. 
# Because read 2 is 35 nt and contains 10 nt sequence to be trimmed anyway. 
# There is no chance for sequencing procedure to incorporate corresponding adapter sequence.
OUTPUT_DATA_REPO = "cutadapt_module"
JULIA_NAME = "ilia_cutadapt_trimmed_reads_11.10.2019.jl"
NAMES = ["read1, read2"]


# Inputs
read_1 = sys.argv[1]  # Command line input 1 for read 1. Assuming it is fastq.gz
read_2 = sys.argv[2]  # Command line input 2 for read 2. Assuming it is fastq.gz
output_dir = sys.argv[3]  # Command line input 3 for the directory for output files, as described in main.py
read_paths = [read_1, read_2]
adapters = [adapter_seq_read_1, adapter_seq_read_2]


# Create and set working directory
julia_script_dir = os.path.join(os.path.dirname(__file__), JULIA_NAME)  # Julia script should be next to this script.
output_dir_module = os.path.join(output_dir, OUTPUT_DATA_REPO)
if not os.access(output_dir_module, os.W_OK) or not os.path.isdir(output_dir_module):  # Create directory if not exist
    os.mkdir(output_dir_module)
os.chdir(output_dir_module)  # Since everything will be output there


# The structure of the final library is as follows:
# **[Read1 sequencing primer annealing seq] - NN - footprint - NNNNN - BBBBB - [Read2 sequencing primer annealing seq]**


# Cutadapt run for adapter trimming
temp_paths = list()
for read_adapter, read_path, read_name in zip(adapters, read_paths, NAMES):  # Run for each read.
    if read_adapter:  # If adapter is defined for a given read. If none, skip it.
        output_for_given_read = f"-o {read_name}_cutadapt_temp.fastq.gz "
        subprocess.run((
            f"{which('cutadapt')} "  # Define which cutadapt installation to use
            f"--cores={cpu_count()} "  # Define how many core to be used. All cores are now using
            f"-q20 -m23 "  # Settings for quality and minimum length to cut the adapter sequence
            "--discard-untrimmed "  # todo:? 
            "-O6 "  # Require minimum length overlap between read and adapter for an adapter to be found.
            f"-a {read_adapter} "
            f"-o {output_for_given_read} "  # Path to output trimmed sequences
            f"{read_path} "  # Input file
            f"1> 'report_cutadapt_temp_{read_name}.txt'"
        ), shell=True)
        temp_paths.append(os.path.join(output_dir_module, output_for_given_read))  # Save the temp output path
    else:
        temp_paths.append(read_path)  # Append the untrimmed read to continue analysis in a tidy way


# Cutadapt run to trim 2 nt UMI at 5' end and 5 nt L1 barcode at 3' end
subprocess.run((
    f"{which('cutadapt')} "  # Define which cutadapt installation to use
    f"--cores={cpu_count()} "  # Define how many core to be used. All cores are now using
    f"-u2 -U5 "  # Settings for trimming 5' and 3' ends
    f"-o {NAMES[0]}_cutadapt.fastq.gz "  # Path to output trimmed sequences for read 1.
    f"-p {NAMES[1]}_cutadapt_untrimmed_umi.fastq.gz "  # Path to output trimmed sequences for read 2. 
    f"{temp_paths[0]} {temp_paths[1]} "  # Input file for 
    f"1> 'report_cutadapt_trimming.txt'"
), shell=True)


# Run Ilia's Julia script to create UMI aware trimmed data
# It should be run only read 2 as it contains the UMI at 5' end
subprocess.run((
    f"{which('julia')} "  # Define which julia installation to use
    f"{julia_script_dir} "  # The path for the script
    f"{NAMES[1]}_cutadapt_untrimmed_umi.fastq.gz "  # Input
    f"{NAMES[1]}_cutadapt_umi_aware.fastq.gz "  # Output
    "--umi5 5"  # Settings
), shell=True)


# End of the script
