"""
Cutadapt for paired end sequencing, runs for both reads
"""


import os
import sys
import joblib
import subprocess
from multiprocessing import cpu_count

repository_name = "RiboSeqProcess"
sys.path.append(os.path.abspath(__file__).split(repository_name)[0] + repository_name)
from archieve.common_functions import *


# Check if necessary packages were installed.
check_exist_package("cutadapt")
check_exist_package("umi_tools")


# SETTINGS
adapter_seq_read_1 = "ATCGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
PATTERN_5prime = "^(?P<umi_1>.{2})"  # Regex to get UMI from first 2 in read1
PATTERN_3prime = "(?P<umi_2>.{5})(?P<discard_2>.{5})$"  # Regex to get rid of barcode and take 5 nt UMI from read 2
PATTERN = PATTERN_5prime + ".*" + PATTERN_3prime


# Inputs
read_1 = sys.argv[1]  # Command line input 1 for read 1. Assuming it is fastq.gz
output_dir = sys.argv[2]  # Command line input 3 for the directory for output files, as described in main.py
temp_dir = sys.argv[3]
OUTPUT_DATA_REPO = "01_preprocessing"


# Create and set working directory
output_dir_module = create_dir(output_dir, OUTPUT_DATA_REPO)  # Create directory if not exist
os.chdir(output_dir_module)  # Since everything will be output there


# The structure of the final library is as follows:
# **[Read1 sequencing primer annealing seq] - NN - footprint - NNNNN - BBBBB - [Read2 sequencing primer annealing seq]**


# Cutadapt run for adapter trimming
temp_path = os.path.join(output_dir_module, "read_1_cutadapt_temp.fastq.gz")  # Outputs for the first run


print("Cutadapt subprocess to remove the adapter sequences is now running.")
subprocess.run((
    f"{which('cutadapt')} "  # Define which cutadapt installation to use
    f"--cores={cpu_count()} "  # Define how many core to be used. All cores are now using
    "--match-read-wildcards "
    f"-q20 -m23 "  # Settings for quality and minimum length to cut the adapter sequence
    # If we don't below two, UMI tool definitely malfunction
    "--discard-untrimmed "
    "-O6 "  # Require minimum length overlap between read and adapter for an adapter to be found.
    f"-a {adapter_seq_read_1} "  # The adapter flags. No flanking white space allowed
    f"-o {temp_path} "  # Path to output trimmed sequences
    f"{read_1} "  # Input file
    f"1> 'report_cutadapt_temp.txt'"
), shell=True)


# Umi-tools
final_path = os.path.join(output_dir_module, "read_1_no-adapt_umi-aware.fastq.gz")

print("Umitool subprocess to remove barcode and extract UMI is now running.")
subprocess.run((
    f"{which('umi_tools')} extract "  # Define which extract installation to use
    "--extract-method=regex "
    f"--bc-pattern='{PATTERN}' "  # Barcode pattern
    f"-I {temp_path} -S {final_path} "  # Input and output for read 1
    f"--log=umi_tools.log" # Log the results
), shell=True)


# Output the final file path to use in a pipeline
joblib.dump(final_path, os.path.join(output_dir, ".01_preprocessing.joblib"))


# End of the script
