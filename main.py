"""
Analysis of Paired End Sequencing of Disome Profiling

The working directory is the directory where the code is run.
The pipeline works only for a pair of reads.
If there is more than o sample or replicates, you have to run the script multiple times.

The output file is created at the working directory.
"""


import os
import subprocess
import sys
import shutil
import multiprocessing


# CONSTANTS
DATA_REPO = "output"  # Output files like sam or fasta
TEMP_REPO = "temp"  # Created files to make scripts work like indexes.


# Inputs
read_1 = sys.argv[1]  # Command line input 1 for read 1
read_2 = sys.argv[2]  # Command line input 2 for read 2


# Create and set working directory
running_directory = os.getcwd() 
output_dir = os.path.join(running_directory, DATA_REPO)
if not os.access(output_dir, os.W_OK) or not os.path.isdir(output_dir):  # Create directory if not exist
    os.mkdir(output_dir)
temp_dir = os.path.join(running_directory, TEMP_REPO)
if not os.access(temp_dir, os.W_OK) or not os.path.isdir(temp_dir):  # Create directory if not exist
    os.mkdir(temp_dir)


# Run Cutadapt module
subprocess.run((
    f"{which('python3')} "  # Define which python installation to use
    "cutadapt.py "  # The script is at the same directory with the main.py
    f"{read_1} "  # sys.argv[1]
    f"{read_2} "  # sys.argv[2]
    f"{output_dir}"  # sys.argv[3]
), shell=True)


# Run rrna_removal.py


# Run link_reads.py


# Run star_alignment.py


