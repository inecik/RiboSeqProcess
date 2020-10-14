"""
The script downloads and filters transcriptome data for Bowtie2 alignment tool.
It filters out everything other than protein coding transcripts.
Then, it runs indexing tool of Bowtie2.
"""


import os
import subprocess
import gzip
from Bio import SeqIO
from datetime import datetime
import re
from multiprocessing import cpu_count


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# Operations for working environment
running_directory = os.getcwd() 
data_repository = "bowtie2_transcriptome"  # name of the database containing folder
data_repo_dir = os.path.join(running_directory, data_repository)
if not os.access(data_repo_dir, os.W_OK) or not os.path.isdir(data_repo_dir): # Create directory if not exist
    print("Data directory created")
    os.mkdir(data_repo_dir)


# Download the transcriptome data
print("Download the transcriptome data")
# Release-96 (GRCh38.p12) is used to be compatible with Mati and Kai's previous works
database_name_fasta = "Homo_sapiens.GRCh38.cdna.all.fa.gz"
subprocess.run((
    f"cd {data_repo_dir};"  # Change the directory to the downloading directory
    "curl -L -R -O ftp://ftp.ensembl.org/pub/release-96/"  # Curl FTP fetch
    f"fasta/homo_sapiens/cdna/{database_name_fasta}"  #.. remaining of the URL
), shell=True)
path_fasta = os.path.join(data_repo_dir, database_name_fasta) # Find out the paths


# Determine the names for the output files
print("Determine the names for the output files")
path_fasta_output = re.search(r".*[^.fa.gz$]", path_fasta).group() + "_filtered.fa"
path_report_output = os.path.join(data_repo_dir, "report" + datetime.now().strftime("-%Y.%m.%d-%H.%M.%S") + ".txt")


# Filter the fasta file
print("Filter the fasta file")
counter_total_fasta = 0  # Count the line in the file
counter_filtered_fasta = 0
with gzip.open(path_fasta, "rt") as handle:  # Open the gz file
    with open(path_fasta_output, "wt") as output_handle:  # Open a file for output
        for record in SeqIO.parse(handle, "fasta"):  # Read the fasta record by a function in Biopython package
            counter_total_fasta += 1
            # Keep only if the transcript is protein coding.
            # Check: https://www.ensembl.org/info/genome/genebuild/biotypes.html
            if  re.search(r"\sgene_biotype:protein_coding\s", record.description):
                SeqIO.write(record, output_handle, "fasta")
                counter_filtered_fasta += 1


# Write down the report
print("Write down the report for filtering:")
fasta_string = f"Original fasta number of line:\t{counter_total_fasta}\n" \
               f"Filtered fasta number of line:\t{counter_filtered_fasta}\n"
with open(path_report_output, "w") as output_handle:
    output_handle.write(f"Filtering Report:\n\n{fasta_string}")


# Run Bowtie2_build function in Bowtie2 module.
print("Run Bowtie2_build function in Bowtie2 module.")
n_core = cpu_count()  # Get the number of cores of the system for multiprocessing
bowtie2_index = os.path.join(data_repository, "bowtie_index")
if not os.access(bowtie2_index, os.W_OK) or not os.path.isdir(bowtie2_index): # Create directory if not exist
    os.mkdir(bowtie2_index)
subprocess.run((
    f"cd {bowtie2_index};"    # Change the directory to the index directory
    "bowtie2-build "  # Name of the function
    f"--threads {n_core} "  # Number of threads to be used
    f"{path_fasta_output} "  # Input file. -f is to indicate the file is in fasta format
    f"{datetime.now().strftime('%Y.%m.%d_%H.%M.%S')}"  # The basename of the index files to write
), shell=True)


#End of the script
