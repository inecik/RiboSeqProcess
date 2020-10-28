"""
It is to remove rRNA, tRNAs among the sequencing reads.
It works with trimmed reads by previous module.
It should accept two pairs as sys.argv
Bowtie2-build should have been done before running this script.
The output folder will be created to the directory where the script is called.
It looks bowtie2_noncoding folder in the directory where the script is called.
"""


import os
import sys
import subprocess
from multiprocessing import cpu_count
from shutil import which


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# CONSTANTS
OUTPUT_DATA_REPO = "bowtie2_rnaremove"  # Name of the database containing folder
INDEX_DIR = "bowtie2_rrna-trna/bowtie2_index/"  # Check database_transcriptome_bowtie2.py
INDEX_BASE = "sapiens-coli_rrna-trna"  # Check database_rnaremove_bowtie2.py
OUTPUT_FASTQ = "norRNA.fastq"
REPORT_FILE = "report_rnaremove.txt"


# Operations for working environment and file name related operations
to_rRNA_depleted_fasta_read1 = sys.argv[1]  # Accept as the first command line argument
to_rRNA_depleted_fasta_read2 = sys.argv[2]  # Accept as the second command line argument
output_dir = sys.argv[3]
temp_dir = sys.argv[4]

data_repo_dir = os.path.join(output_dir, OUTPUT_DATA_REPO)
if not os.access(data_repo_dir, os.W_OK) or not os.path.isdir(data_repo_dir):  # Create directory if not exist
    os.mkdir(data_repo_dir)

index_files = os.path.join(temp_dir, INDEX_DIR, INDEX_BASE)
output_report_path = os.path.join(data_repo_dir, REPORT_FILE)
output_path = os.path.join(data_repo_dir, OUTPUT_FASTQ)

bowtie2_run = subprocess.run((
    f"cd {data_repo_dir};"  # Change the directory to the index directory
    f"{which('bowtie2')} "  # Run Bowtie2 module
    f"-p{cpu_count()} "  # Number of core to use
    "-D30 -R5 -L13 " # D and R added by me. L is from the manuscript.
    "-q "  # Indicates the inputs are fastq
    "--no-mixed "  # Do not search for individual pairs if one in a pair does not align.
    "--time "  # Print the wall-clock time required to load the index files and align the reads. 
    f"-x {index_files} "  # Index directory with the base name
    f"-1 {to_rRNA_depleted_fasta_read1} "  # Read 1
    f"-2 {to_rRNA_depleted_fasta_read2} "  # Read 2
    f"--un-conc {output_path} "  # Output fastq file, Contains all reads which did not aligned RNAs.
    "-S /dev/null "  # Discard alignment sam file
    f"2> {output_report_path}"
), shell=True)


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
