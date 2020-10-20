"""
Cutadapt for paired end sequencing, runs for both reads
"""


import os
import sys
import subprocess
from multiprocessing import cpu_count
from shutil import which


# SETTINGS
adapter_seq_read_1 = "ATCGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
adapter_seq_read_2 = None
# Mati and Kai's protocol does not necessitates to cut read 2.
# Because read 2 is 35 nt and contains 10 nt sequence to be trimmed anyway.
# There is no chance for sequencing procedure to incorporate corresponding adapter sequence.


# Inputs
read_1 = sys.argv[1]  # Command line input 1 for read 1. Assuming it is fastq.gz
read_2 = sys.argv[2]  # Command line input 2 for read 2. Assuming it is fastq.gz
output_dir = sys.argv[3]  # Command line input 3 for the directory for output files, as described in main.py
read_paths = [read_1, read_2]
NAMES = ["read1", "read2"]
OUTPUT_DATA_REPO = "cutadapt_module"
JULIA_NAME = "ilia_cutadapt_trimmed_reads_11.10.2019.jl"
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
temp_paths = [f"{NAMES[0]}_cutadapt_temp.fastq.gz", f"{NAMES[1]}_cutadapt_temp.fastq.gz"]  # Outputs for the first run
temp_paths = [os.path.join(output_dir_module, i) for i in temp_paths]
read1_adapter = f"-a {adapters[0]}" if adapters[0] else ""  # Create flag is adapter is provided
read2_adapter = f"-A {adapters[1]}" if adapters[1] else ""  # Create flag is adapter is provided

subprocess.run((
    f"{which('cutadapt')} "  # Define which cutadapt installation to use
    f"--cores={cpu_count()} "  # Define how many core to be used. All cores are now using
    f"-q20 -m23 "  # Settings for quality and minimum length to cut the adapter sequence
    # "--discard-untrimmed "  # todo: Understand/Ask why!?!
    # "-O6 "  # Require minimum length overlap between read and adapter for an adapter to be found.  # todo: Understand/Ask why!?!
    f"{read1_adapter} {read2_adapter}".strip() + " "  # The adapter flags. No flanking white space allowed
    f"-o {temp_paths[0]} -p {temp_paths[1]} "  # Path to output trimmed sequences
    f"{read_paths[0]} {read_paths[1]} "  # Input file
    f"1> 'report_cutadapt_temp.txt'"
), shell=True)


# Cutadapt run to trim 2 nt UMI at 5' end and 5 nt L1 barcode at 3' end
temp2_paths = [f"{NAMES[0]}_cutadapt.fastq.gz", f"{NAMES[1]}_cutadapt_untrimmed_umi.fastq.gz"]
temp2_paths = [os.path.join(output_dir_module, i) for i in temp2_paths]

subprocess.run((
    f"{which('cutadapt')} "  # Define which cutadapt installation to use
    f"--cores={cpu_count()} "  # Define how many core to be used. All cores are now using
    f"-u2 -U5 "  # Settings for trimming 5' and 3' ends
    f"-o {temp2_paths[0]} -p {temp2_paths[1]} "  # Path to output trimmed sequences
    f"{temp_paths[0]} {temp_paths[1]} "  # Input file for 
    f"1> 'report_cutadapt_trimming.txt'"
), shell=True)


# Run Ilia's Julia script to create UMI aware trimmed data
# It should be run only read 2 as it contains the UMI at 5' end
umi_output_path = os.path.join(output_dir_module, f"{NAMES[1]}_cutadapt_umi_aware.fastq.gz")
subprocess.run((
    f"{which('julia')} "  # Define which julia installation to use
    f"{julia_script_dir} "  # The path for the script
    f"{temp2_paths[1]} "  # Input
    f"{umi_output_path} "  # Output
    "--umi5 5"  # Settings
), shell=True)


# End of the script
