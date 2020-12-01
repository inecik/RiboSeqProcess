"""
Star alignment for single end sequencing.

"""


import os
import sys
import subprocess
from multiprocessing import cpu_count
from shutil import which

repository_name = "RiboSeqProcess"
sys.path.append(os.path.abspath(__file__).split(repository_name)[0] + repository_name)
from archieve.common_functions import *


# Check if necessary packages were installed.
check_exist_package("STAR")


# CONSTANTS
OUTPUT_DATA_REPO = "04_genomealignment"
TEMP_DATA_REPO = "04_genomealignment"
GENOME_INDEX_DIR_NAME = "index_genome"


# Inputs
read_file = sys.argv[1]  # For read 1
output_dir = sys.argv[2]  # For the directory for output files, as described in main.py
temp_dir = sys.argv[3]  # For the directory for temp files, as described in main.py


# Create and set working directory for genome index creation
temp_dir_module = os.path.join(temp_dir, TEMP_DATA_REPO)
genome_index_dir = os.path.join(temp_dir_module, GENOME_INDEX_DIR_NAME)


# Create and set working directory for alignment
output_dir_module = create_dir(output_dir, OUTPUT_DATA_REPO)
os.chdir(output_dir_module)  # Since everything will be output there


# Actual alignment for single end sequencing
subprocess.run((
    f"{which('STAR')} "  # Define which star installation to use
    f"--runThreadN {cpu_count()} "  # Define how many core to be used. All cores are now using
    f"--genomeDir {genome_index_dir} "  # Directory for genome index which has been just created
    f"--readFilesIn {read_file} "
    # All parameters were inherited from Mati-Kai's pipeline.
    "--outFilterMultimapNmax 1 "  
    # "--peOverlapNbasesMin 6 " It is for paired end sequencing
    # "--peOverlapMMp 0.1 " It is for paired end sequencing
    "--outFilterType BySJout "
    "--alignIntronMin 5 "
    f"--outFileNamePrefix {os.path.join(output_dir_module, 'single_')} "
    "--outReadsUnmapped Fastx "
    "--outSAMtype BAM SortedByCoordinate "
    "--outSAMattributes All XS "
    "--quantMode GeneCounts "
    "--twopassMode Basic"
), shell=True)


# End of the script
