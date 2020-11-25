"""
Star alignment for paired end sequencing.

"""


import os
import sys
import subprocess
from multiprocessing import cpu_count
from shutil import which


# CONSTANTS
OUTPUT_DATA_REPO = "star_paired_module"
TEMP_DATA_REPO = "star_module"
GENOME_INDEX_DIR_NAME = "genome_index"


# Inputs
read_1 = sys.argv[1]  # Command line input 1 for read 1. Assuming it is fastq.gz
read_2 = sys.argv[2]  # Command line input 2 for read 2. Assuming it is fastq.gz
output_dir = sys.argv[3]  # Command line input 3 for the directory for output files, as described in main.py
temp_dir = sys.argv[4]  # Command line input 4 for the directory for temp files, as described in main.py


# Create and set working directory for genome index creation
temp_dir_module = os.path.join(temp_dir, TEMP_DATA_REPO)
genome_index_dir = os.path.join(temp_dir_module, GENOME_INDEX_DIR_NAME)


# Create and set working directory for alignment
output_dir_module = os.path.join(output_dir, OUTPUT_DATA_REPO)
if not os.access(output_dir_module, os.W_OK) or not os.path.isdir(output_dir_module):  # Create directory if not exist
    os.mkdir(output_dir_module)
os.chdir(output_dir_module)  # Since everything will be output there


# Actual alignment for paired end sequencing
subprocess.run((
    f"{which('STAR')} "  # Define which star installation to use
    f"--runThreadN {cpu_count()} "  # Define how many core to be used. All cores are now using
    f"--genomeDir {genome_index_dir} "  # Directory for genome index which has been just created
    f"--readFilesIn {read_1} {read_2} "
    # All parameters were inherited from Mati-Kai's pipeline.
    "--outFilterMultimapNmax 1 "  
    "--peOverlapNbasesMin 6 "
    "--peOverlapMMp 0.1 "
    "--outFilterType BySJout "
    "--alignIntronMin 5 "
    f"--outFileNamePrefix {os.path.join(output_dir_module, 'paired_')} "
    "--outReadsUnmapped Fastx "
    "--outSAMtype BAM SortedByCoordinate "
    "--outSAMattributes All XS "
    "--quantMode GeneCounts "
    "--twopassMode Basic"
), shell=True)


# End of the script
