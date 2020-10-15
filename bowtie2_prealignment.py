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
running_directory = os.getcwd()
data_repo_dir = os.path.join(running_directory, DATA_REPO)
if not os.access(data_repo_dir, os.W_OK) or not os.path.isdir(data_repo_dir):  # Create directory if not exist
    print("Data directory created")
    os.mkdir(data_repo_dir)
output_path = os.path.join(data_repo_dir, SAM_FILE)
index_files = os.path.join(running_directory, INDEX_DIR, INDEX_BASE)
rRNA_depleted_fasta_read1 = sys.argv[1]  # Accept as the first command line argument
rRNA_depleted_fasta_read2 = sys.argv[2]  # Accept as the second command line argument


bowtie2_run = subprocess.run((
    f"cd {data_repo_dir};"  # Change the directory to the index directory
    f"{which('bowtie2')} "  # Run Bowtie2 module
    "-D20 -R3"  # Increases sensitivity 
    "-I35 -X85 "  # Search only those that has 35-85 nt. Makes Bowtie2 slower. 
    # Todo: Increase 85 to 100 and decrease 35 to 20. Rerun everthing and check the footprints
    "-q "  # Specifies the inputs are in fastq format 
    f"-p{cpu_count()} "  # Number of core to use
    "--no-discordant "  # Filter pairs does not obey orientation/length constraints 
    "--no-mixed "  # Do not search for individual pairs if one in a pair does not align.
    "--end-to-end "  # Indicate the alignment is end to end. No soft clipping will be applied.
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
