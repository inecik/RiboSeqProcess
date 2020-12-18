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
import re
from multiprocessing import cpu_count
from shutil import which
from Bio import SeqIO
import joblib

repository_name = "RiboSeqProcess"
sys.path.append(os.path.abspath(__file__).split(repository_name)[0] + repository_name)
from archieve.common_functions import *


# Check if necessary packages were installed.
check_exist_package("bowtie2-build")
ensembl_release = 102


# CONSTANTS
DATA_REPO = "03_linkpairs"  # name of the database containing folder
DATABASE_FASTA = "Homo_sapiens.GRCh38.cdna.all.fa.gz"  # Name of fasta file in the ftp server
INDEX_BASE = "homo_sapiens_protein_coding_transcriptome"  # The basename of the index files to write

# Operations for working environment
temp_dir = sys.argv[1]
temp_module_dir = create_dir(temp_dir, DATA_REPO)


# Download the transcriptome data and file name related operations
path_fasta = os.path.join(temp_module_dir, DATABASE_FASTA)  # Find out the paths
if not os.access(path_fasta, os.R_OK) or not os.path.isfile(path_fasta):
    print("Database for human transcriptome is now downloading.")
    subprocess.run((
        f"cd {temp_module_dir};"  # Change the directory to the downloading directory
        f"curl -L -R -O ftp://ftp.ensembl.org/pub/release-{ensembl_release}/"  # Curl FTP fetch
        f"fasta/homo_sapiens/cdna/{DATABASE_FASTA}; "
        f"gzip -d {path_fasta}"), shell=True)  # .. remaining of the URL
path_fasta_output = re.search(r"(.*)\.fa\.gz$", path_fasta).group(1) + ".fa"


# Run Bowtie2_build function in Bowtie2 module.
bowtie2_index = create_dir(temp_module_dir, "index_transcriptome")
metadata_index_path = os.path.join(temp_module_dir, ".03_linkpairs_metadata.joblib")

try:
    # Check if the index file
    assert os.path.isfile(metadata_index_path) and os.access(metadata_index_path, os.R_OK)
    metadata_index_previously = joblib.load(metadata_index_path)
    metadata_index = get_files_metadata(bowtie2_index)
    assert metadata_index_previously == metadata_index
    print("Index files for link pairing are already present.")

except:
    print("Index files for link pairing are now being created.")
    subprocess.run((
        f"cd {bowtie2_index}; "  # Change the directory to the index directory
        f"{which('bowtie2-build')} "  # Name of the function
        f"--threads {cpu_count()} "  # Number of threads to be used.
        f"{path_fasta_output} "  # Input file. -f is to indicate the file is in fasta format
        f"{INDEX_BASE} "  # The basename of the index files to write
        "> report_bowtie2_database.log"
    ), shell=True)

    # Write the file info
    metadata_index = get_files_metadata(bowtie2_index)
    joblib.dump(metadata_index, metadata_index_path)


# End of the script
