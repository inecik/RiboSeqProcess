"""
The script downloads and filters transcriptome data for Bowtie2 alignment tool.
It filters out everything other than protein coding transcripts.
Then, it runs indexing tool of Bowtie2.
The output folder will be created to the directory where the script is called.
"""


import gzip
import os
import re
import subprocess
import sys
from datetime import datetime
from multiprocessing import cpu_count
from shutil import which
from Bio import SeqIO

sys.path.append(os.path.abspath(__file__).split("Ribo-seq-Analysis")[0] + "Ribo-seq-Analysis")
from supplementary.common_functions import bcolors as c


# Check if necessary packages were installed.
if not which('bowtie2-build'):
    sys.exit(f"{c.FAIL}Bowtie2 package should be installed.{c.ENDC}")
else:
    print("Bowtie2 installation is found.")


# CONSTANTS
DATA_REPO = "bowtie2_transcriptome"  # name of the database containing folder
DATABASE_FASTA = "Homo_sapiens.GRCh38.cdna.all.fa.gz"  # Name of fasta file in the ftp server
INDEX_BASE = "homo_sapiens_protein_coding_transcriptome"  # The basename of the index files to write

# Operations for working environment
temp_dir = sys.argv[1]
temp_module_dir = os.path.join(temp_dir, DATA_REPO)
if not os.access(temp_module_dir, os.W_OK) or not os.path.isdir(temp_module_dir):  # Create directory if not exist
    os.mkdir(temp_module_dir)
    print("Temp directory created.")


# Download the transcriptome data and file name related operations
print("Database for human transcriptome is now downloading.")
subprocess.run((
    f"cd {temp_module_dir};"  # Change the directory to the downloading directory
    "curl -L -R -O ftp://ftp.ensembl.org/pub/release-100/"  # Curl FTP fetch
    f"fasta/homo_sapiens/cdna/{DATABASE_FASTA}"), shell=True)  # .. remaining of the URL
path_fasta = os.path.join(temp_module_dir, DATABASE_FASTA)  # Find out the paths
path_fasta_output = re.search(r"(.*)\.fa\.gz$", path_fasta).group(1) + "_filtered.fa"
path_report_output = os.path.join(temp_module_dir, "report" + datetime.now().strftime("-%Y.%m.%d-%H.%M.%S") + ".txt")


# Filter the fasta file
print("Database for human transcriptome is being filtered for protein coding RNAs.")
counter_total_fasta = 0  # Count the line in the file
counter_filtered_fasta = 0
with gzip.open(path_fasta, "rt") as handle:  # Open the gz file
    with open(path_fasta_output, "wt") as output_handle:  # Open a file for output
        for record in SeqIO.parse(handle, "fasta"):  # Read the fasta record by a function in Biopython package
            counter_total_fasta += 1
            # Keep only if the transcript is protein coding.
            # Check: https://www.ensembl.org/info/genome/genebuild/biotypes.html
            if re.search(r"\stranscript_biotype:protein_coding\s", record.description):
                SeqIO.write(record, output_handle, "fasta")
                counter_filtered_fasta += 1


# Write down the report
fasta_string = f"Original fasta number of line:\t{counter_total_fasta}\n" \
               f"Filtered fasta number of line:\t{counter_filtered_fasta}\n"
with open(path_report_output, "w") as output_handle:
    output_handle.write(f"Filtering Report:\n\n{fasta_string}")


print("Human transcriptome indexing is started.")
# Run Bowtie2_build function in Bowtie2 module.
n_core = cpu_count()  # Get the number of cores of the system for multiprocessing
bowtie2_index = os.path.join(temp_module_dir, "bowtie2_index")
if not os.access(bowtie2_index, os.W_OK) or not os.path.isdir(bowtie2_index):  # Create directory if not exist
    os.mkdir(bowtie2_index)
subprocess.run((
    f"cd {bowtie2_index}; "  # Change the directory to the index directory
    f"{which('bowtie2-build')} "  # Name of the function
    f"--threads {n_core} "  # Number of threads to be used.
    f"{path_fasta_output} "  # Input file. -f is to indicate the file is in fasta format
    f"{INDEX_BASE} "  # The basename of the index files to write
    "> report_bowtie2_database.log"
), shell=True)


# End of the script
