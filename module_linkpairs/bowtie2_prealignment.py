"""
It is to align rRNA depleted and trimmed reads to the transcriptome data.
Accepts read 1 and read 2 as sys.argv respectively.
Bowtie2-build should have been done before running this script.
The output folder will be created to the directory where the script is called.
It looks bowtie2_transcriptome folder in the directory where the script is called.
"""


import os
import sys
import subprocess
from multiprocessing import cpu_count
from datetime import datetime
from shutil import which


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# CONSTANTS
DATA_REPO = "bowtie2_prealignment"  # Name of the database containing folder
SAM_FILE = "bowtie2_prealignment"
INDEX_DIR = "bowtie2_transcriptome/bowtie2_index/"  # Check database_transcriptome_bowtie2.py
INDEX_BASE = "homo_sapiens_protein_coding_transcriptome"  # Check database_transcriptome_bowtie2.py


# Operations for working environment and file name related operations
rRNA_depleted_fasta_read1 = sys.argv[1]  # Accept as the first command line argument
rRNA_depleted_fasta_read2 = sys.argv[2]  # Accept as the second command line argument
output_dir = sys.argv[3]
temp_dir = sys.argv[4]
output_module_dir = os.path.join(output_dir, DATA_REPO)
if not os.access(output_module_dir, os.W_OK) or not os.path.isdir(output_module_dir):  # Create directory if not exist
    os.mkdir(output_module_dir)
output_path = os.path.join(output_module_dir, SAM_FILE)
index_files = os.path.join(temp_dir, INDEX_DIR, INDEX_BASE)


bowtie2_run = subprocess.run((
    f"cd {output_module_dir};"  # Change the directory to the index directory
    f"{which('bowtie2')} "  # Run Bowtie2 module
    "-D 40 -R 6 -N 0 -L 20 -i S,1,0.50 "  # Alignment effort and sensitivity. It is now very high.
    "-I20 -X150 "  # Search only those that has 30-180 nt. Makes Bowtie2 slower. 
    "--score-min G,20,6 "  # min score lowered
    "--ma 3 "  # ma bonus increased
    "--no-unal "  # To suppress the non-aligned reads
    "-q "  # Specifies the inputs are in fastq format 
    f"-p{cpu_count()} "  # Number of core to use
    "--no-discordant "  # Filter pairs does not obey orientation/length constraints 
    "--no-mixed "  # Do not search for individual pairs if one in a pair does not align.
    "--local "  # Indicate the alignment is local. Soft clipping will be applied.
    "--time "  # Print the wall-clock time required to load the index files and align the reads. 
    f"-x {index_files} "  # Index directory with the base name
    f"-1 {rRNA_depleted_fasta_read1} "  # Read 1
    f"-2 {rRNA_depleted_fasta_read2} "  # Read 2
    f"-S {output_path}.sam"  # Output sam file
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
