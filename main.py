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


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# CONSTANTS
OUTP_MAIN = "OUTPUT"  # Output files like sam or fasta
TEMP_MAIN = "TEMP"  # Created files to make scripts work like indexes.
FIGU_MAIN = "FIGURES"


# Inputs
raw_inputs = {'read_1': sys.argv[1], 'read_2': sys.argv[2]}  # Command line inputs for the reads


# Obtain the root dir names for scripts to run properly
running_directory = os.getcwd()  # Where this script is called
scripts_directory = os.path.dirname(__file__)  # Where this package is


# Create directories to work on
output_dir = os.path.join(running_directory, OUTP_MAIN)
if not os.access(output_dir, os.W_OK) or not os.path.isdir(output_dir):  # Create directory if not exist
    os.mkdir(output_dir)

tempor_dir = os.path.join(running_directory, TEMP_MAIN)
if not os.access(tempor_dir, os.W_OK) or not os.path.isdir(tempor_dir):  # Create directory if not exist
    os.mkdir(tempor_dir)

figure_dir = os.path.join(running_directory, FIGU_MAIN)
if not os.access(figure_dir, os.W_OK) or not os.path.isdir(figure_dir):  # Create directory if not exist
    os.mkdir(figure_dir)


# ____________________________________________________________
# Cutadapt Module
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾


subprocess.run((
    f"{which('python3')} "  # Define which python installation to use
    f"{os.path.join(scripts_directory, "module_cutadapt/cutadapt.py")} "  # Script directory relative to main.py
    f"{raw_inputs['read_1']} "  # sys.argv[1] 
    f"{raw_inputs['read_2']} "  # sys.argv[2]
    f"{output_dir}"  # sys.argv[3]
), shell=True)


# ____________________________________________________________ 
# RNA Remove Module
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

rna_remove_inputs = {'read_1': "read1_cutadapt.fastq.gz", 'read_2': "read2_cutadapt_umi_aware.fastq.gz"}

# Genome indexes
subprocess.run((  # todo: check if exist, run otherwise
    f"{which('python3')} "  # Define which python installation to use
    f"{os.path.join(scripts_directory, "module_rnaremove/database_rnaremove_bowtie2.py")} "  # Script directory relative to main.py
    f"{TEMP_MAIN} "  # sys.argv[1] 
), shell=True)

# Alignment

#### todo INCOMPLETE HERE
subprocess.run((
    f"{which('python3')} "  # Define which python installation to use
    f"{os.path.join(scripts_directory, "module_rnaremove/bowtie2_rnaremove.py")} "  # Script directory relative to main.py
    f"{rna_remove_inputs['read_1']} "  # sys.argv[1] 
    f"{rna_remove_inputs['read_2']} "  # sys.argv[2]
    f"{output_dir}"  # sys.argv[3]
), shell=True)

# Run link_reads.py


# Run star_alignment.py


# Run assignment