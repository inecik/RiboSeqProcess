"""
The script downloads and filters transcriptome data.
It prepares the database for alignment tools.
"""


import os
import subprocess
import gzip
from Bio import SeqIO
from datetime import datetime
import re


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# Operations for working environment
running_directory = os.getcwd() 
data_repository = "data"  # name of the database containing folder
data_repo_dir = os.path.join(running_directory, data_repository)
if not os.access(running_directory, os.W_OK) or not os.path.isdir(data_repo_dir): # Create directory if not exist
    print("Data directory created")
    os.mkdir(data_repo_dir)


# Download the transcriptome data
print("Download the transcriptome data")
# Release-96 (GRCh38.p12) is used to be compatible with Mati and Kai's previous works
database_name_fasta = "Homo_sapiens.GRCh38.cdna.all.fa.gz"
database_name_gtf = "Homo_sapiens.GRCh38.96.chr_patch_hapl_scaff.gtf.gz"
bash_command_download_fasta = (f"cd {data_repo_dir}; curl -L -R -O ftp://ftp.ensembl.org/pub/release-96/"
                                  f"fasta/homo_sapiens/cdna/{database_name_fasta}")
bash_command_download_gtf = (f"cd {data_repo_dir}; curl -L -R -O ftp://ftp.ensembl.org/pub/release-96/"
                                f"gtf/homo_sapiens/{database_name_gtf}")
# Actual download starts
subprocess.run(bash_command_download_fasta, shell=True)
subprocess.run(bash_command_download_gtf, shell=True)
# Find out the paths
path_fasta = os.path.join(data_repo_dir, database_name_fasta)
path_gtf = os.path.join(data_repo_dir, database_name_gtf)


# Determine the names for the output files
print("Determine the names for the output files")
path_fasta_output = re.search(r".*[^.fa.gz$]", path_fasta).group() + "_filtered.fa"
path_gtf_output = re.search(r".*[^.gtf.gz$]", path_gtf).group() + "_filtered.gtf"
path_report_output = os.path.join(data_repo_dir, "report" + datetime.now().strftime("-%Y.%m.%d-%H.%M.%S") + ".txt")


# Filter the fasta file
print("Filter the fasta file")
longest_isoform = dict()  # Keep track of longest transcripts in a given gene
counter_total_fasta = 0  # Count the line in the file
with gzip.open(path_fasta, "rt") as handle:  # Open the gz file
    for record in SeqIO.parse(handle, "fasta"):  # Read the fasta record by a function in Biopython package
        # Keep only if the transcript is protein coding.
        # Check: https://www.ensembl.org/info/genome/genebuild/biotypes.html
        if  re.search(r"\sgene_biotype:protein_coding\s", record.description):
            gene = re.search(r"gene:ENSG[0-9.]*", record.description).group()  # Fetch the associated gene by regex
            counter_total_fasta += 1
            if gene not in longest_isoform:  # if there is not gene in the dict, keep the new record
                longest_isoform[gene] = record
            elif len(longest_isoform[gene].seq) < len(record.seq):  # if the recorded transcript is shorter, update
                longest_isoform[gene] = record


# Write the filtered results to compressed fasta file
print("Write the filtered results to compressed fasta file")
with open(path_fasta_output, "wt") as output_handle:
    SeqIO.write(longest_isoform.values(), output_handle, "fasta")


# Apply the filter to corresponding gtf file
print("Apply the filter to corresponding gtf file")
kept_transcripts = set() # Create a set with transcript Ensembl IDs.
for record in longest_isoform.values():
    # Regex is to remove the '.' if exist, since GTF file doesn't have it
    kept_transcripts.add(re.search(r"ENST[0-9]*", record.id).group())


counter_total_gtf = 1  # Count the line in the file
counter_new_gtf = 0  # Count the lines passing the filter
with gzip.open(path_gtf, "rt") as handle:  # Open gtf.gz file
    with open(path_gtf_output, "w") as output_handle:
        line = handle.readline()
        while line:  # Read line by line until the end
            counter_total_gtf += 1
            if line[0] == "#":  # Do not touch if it is title line
                output_handle.write(line)
                counter_new_gtf += 1
            else:
                id = re.search(r"transcript_id \"(ENST[0-9]*)\"", line)  # Get the transcript info of the entry
                if id and id.group(1) in kept_transcripts:  # Keep only the ones that has transcript_id and in the set
                    output_handle.write(line)
                    counter_new_gtf += 1
            line = handle.readline()  # Go to the next line


# Write down the report
print("Write down the report")
fasta_string = f"Original fasta number of line:\t{counter_total_fasta}\n" \
               f"Filtered fasta number of line:\t{len(longest_isoform)}\n"
gtf_string   = f"Original gtf number of line:\t{counter_total_gtf}\n" \
               f"Filtered gtf number of line:\t{counter_new_gtf}\n"
with open(path_report_output, "w") as output_handle:
    output_handle.write(f"REPORT:\n\n{fasta_string}\n{gtf_string}")


#End of the script
