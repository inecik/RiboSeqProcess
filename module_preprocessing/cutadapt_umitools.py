"""
Cutadapt for paired end sequencing, runs for both reads
"""


import os
import sys
import joblib
import subprocess
from multiprocessing import cpu_count
from shutil import which

sys.path.append(os.path.abspath(__file__).split("Ribo-seq-Analysis")[0] + "Ribo-seq-Analysis")
from module_supplementary.common_functions import bcolors as c


# Check if necessary packages were installed.
if not which('cutadapt') or not which('umi_tools'):
    sys.exit(f"{c.FAIL}Cutadapt and/or umi_tools packages should be installed.{c.ENDC}")
else:
    print("Cutadapt and umi_tools installations are found.")


# SETTINGS
adapter_seq_read_1 = "ATCGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
adapter_seq_read_2 = None
# Mati and Kai's protocol does not necessitates to cut read 2.
# Because read 2 is 35 nt and contains 10 nt sequence to be trimmed anyway.
# There is no chance for sequencing procedure to incorporate corresponding adapter sequence.
PATTERN_read1 = "^(?P<umi_1>.{2}).*"  # Regex to get UMI from first 2 in read1
PATTERN_read2 = ".*(?P<umi_2>.{5})(?P<discard_2>.{5})$"  # # Regex to get rid of barcode and take 5 nt UMI from read 2


# Inputs
read_1 = sys.argv[1]  # Command line input 1 for read 1. Assuming it is fastq.gz
read_2 = sys.argv[2]  # Command line input 2 for read 2. Assuming it is fastq.gz
output_dir = sys.argv[3]  # Command line input 3 for the directory for output files, as described in main.py
temp_dir = sys.argv[4]
read_paths = [read_1, read_2]
NAMES = ["read_1", "read_2"]
OUTPUT_DATA_REPO = "preprocessing_module"
adapters = [adapter_seq_read_1, adapter_seq_read_2]


# Create and set working directory
output_dir_module = os.path.join(output_dir, OUTPUT_DATA_REPO)
if not os.access(output_dir_module, os.W_OK) or not os.path.isdir(output_dir_module):  # Create directory if not exist
    os.mkdir(output_dir_module)
    print("Output directory created.")
os.chdir(output_dir_module)  # Since everything will be output there


# The structure of the final library is as follows:
# **[Read1 sequencing primer annealing seq] - NN - footprint - NNNNN - BBBBB - [Read2 sequencing primer annealing seq]**


# Cutadapt run for adapter trimming
temp_paths = [f"{NAMES[0]}_cutadapt_temp.fastq.gz", f"{NAMES[1]}_cutadapt_temp.fastq.gz"]  # Outputs for the first run
temp_paths = [os.path.join(output_dir_module, i) for i in temp_paths]
read1_adapter = f"-a {adapters[0]}" if adapters[0] else ""  # Create flag is adapter is provided
read2_adapter = f"-A {adapters[1]}" if adapters[1] else ""  # Create flag is adapter is provided


print("Cutadapt subprocess to remove the adapter sequences is now running.")
subprocess.run((
    f"{which('cutadapt')} "  # Define which cutadapt installation to use
    f"--cores={cpu_count()} "  # Define how many core to be used. All cores are now using
    "--match-read-wildcards "
    f"-q20 -m23 "  # Settings for quality and minimum length to cut the adapter sequence
    # "--discard-untrimmed "  # todo: Understand/Ask why!?!
    # "-O6 "  # Require minimum length overlap between read and adapter for an adapter to be found.  # todo: Understand/Ask why!?!
    f"{read1_adapter} {read2_adapter}".strip() + " "  # The adapter flags. No flanking white space allowed
    f"-o {temp_paths[0]} -p {temp_paths[1]} "  # Path to output trimmed sequences
    f"{read_paths[0]} {read_paths[1]} "  # Input file
    f"1> 'report_cutadapt_temp.txt'"
), shell=True)


# Umi-tools
final_paths = [f"{i}_no-adapt_umi-aware.fastq.gz" for i in NAMES]
final_paths = [os.path.join(output_dir_module, i) for i in final_paths]

print("Umitool subprocess to remove barcode and extract UMI is now running.")
subprocess.run((
    f"{which('umi_tools')} extract "  # Define which extract installation to use
    "--extract-method=regex "
    f"--bc-pattern='{PATTERN_read1}' "  # Barcode pattern
    f"--bc-pattern2='{PATTERN_read2}' "  # Barcode pattern for paired reads"
    f"-I {temp_paths[0]} -S {final_paths[0]} "  # Input and output for read 1
    f"--read2-in={temp_paths[1]} --read2-out={final_paths[1]} "  # Input and output for read 2
    f"--log=umi_tools.log" # Log the results
), shell=True)


# Output the final file path to use in a pipeline
joblib.dump(final_paths, os.path.join(temp_dir, ".module_preprocessing_paths.joblib"))


# End of the script
