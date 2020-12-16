"""
It is to remove rRNA, tRNAs among the sequencing reads.
It works with trimmed reads by previous module.
It should accept two pairs as sys.argv
Bowtie2-build should have been done before running this script.
If a fifth argument is given 'aligned_fastq' then the aligned reads also will be written to the disk.
"""


import os
import sys
import subprocess
from multiprocessing import cpu_count
from shutil import which
import joblib

repository_name = "RiboSeqProcess"
sys.path.append(os.path.abspath(__file__).split(repository_name)[0] + repository_name)
from archieve.common_functions import *


# Check if necessary packages were installed.
check_exist_package("bowtie2")


# CONSTANTS
OUTPUT_DATA_REPO = "02_cleanup"  # Name of the database containing folder
INDEX_DIR = "02_cleanup/index_ncrna/"  # Check database_transcriptome_bowtie2.py
INDEX_BASE = "homo_sapiens_rrna_trna"  # Check database_rnaremove_bowtie2.py
OUTPUT_FASTQ = "Read_%_norRNA.fastq"



# Operations for working environment and file name related operations
input_read_1 = sys.argv[1]  # Accept as the first command line argument
input_read_2 = sys.argv[2]  # Accept as the second command line argument
output_dir = sys.argv[3]
temp_dir = sys.argv[4]


# Directory settings
data_repo_dir = create_dir(output_dir, OUTPUT_DATA_REPO)
index_files = os.path.join(temp_dir, INDEX_DIR, INDEX_BASE)
output_path = os.path.join(data_repo_dir, OUTPUT_FASTQ)


# rRNA removal
print("Alignment for rRNA removal is now started.")
bowtie2_run = subprocess.run((
    f"cd {data_repo_dir};"  # Change the directory to the index directory
    f"{which('bowtie2')} "  # Run Bowtie2 module
    f"-p{cpu_count()} "  # Number of core to use
    "--no-mixed "  # Do not search for individual pairs if one in a pair does not align.
    "-I20 -X120 " # Default -I=0, -X=500. Since I will disregard X>120 and I<20 in the link-pairing module
    "--time "  # Print the wall-clock time required to load the index files and align the reads.
    "--score-min G,20,6 "  # Allow looser alignment. --ma 3 
    "--local --sensitive-local "  # Apply soft clipping when aligning. Default sensitivity.
    "-k1 "  # We are not interested in the best alignment as long as it aligns somewhere in the indexed fasta.
    f"-x {index_files} "  # Index directory with the base name
    "-q "  # Indicates the inputs are fastq
    f"-1 {input_read_1} "  # Read 1
    f"-2 {input_read_2} "  # Read 2
    f"--un-conc {output_path} "  # Output fastq file, Contains all reads which did not aligned RNAs. 
    "-S /dev/null "  # Discard alignment sam file /dev/null
    f"2> report_rnaremove.txt"
), shell=True)


# Output the final file path to use in a pipeline
final_paths = [os.path.join(data_repo_dir, j) for j in ["Read_1_norRNA.fastq", "Read_2_norRNA.fastq"]]
joblib.dump(final_paths, os.path.join(output_dir, ".02_cleanup.joblib"))


# End of the script

# TROUBLESHOOTING:
#
# Error message:
#
# perl: warning: Setting locale failed.
# perl: warning: Please check that your locale settings:
#     LANGUAGE = (unset),
#     LC_ALL = (unset),
#     LANG = "en_US.UTF-8"
# are supported and installed on your system.
# perl: warning: Falling back to the standard locale ("C").
#
# Solution:
#
# Write the following lines to_bash_profile.
# export LANGUAGE=en_US.UTF-8
# export LC_ALL=en_US.UTF-8
# export LANG=en_US.UTF-8
# export LC_CTYPE=en_US.UTF-8
