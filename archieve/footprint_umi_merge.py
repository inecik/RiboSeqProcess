"""
It is to bring back UMI information to the output of sam_processor.py
It looks bowtie2_prealignment folder in the directory where the script is called.
The output folder will be mentioned bowtie2_prealignment folder.
"""


import os
from Bio import SeqIO
import sys
import subprocess
from multiprocessing import cpu_count
from shutil import which



# CONSTANTS
DATA_REPO = "bowtie2_prealignment"  # Name of the database containing folder
FASTA_FILE = "footprints.fasta"  # Check sam_processor.py
OUTPUT_FILE = "footprints_UMI.fasta"


# Operations for working environment and file name related operations
read_with_umi = sys.argv[1]  # Important: Only FastQ, which contains UMI info in the title
output_dir = sys.argv[2]
output_module_dir = os.path.join(output_dir, DATA_REPO)
input_path = os.path.join(output_module_dir, FASTA_FILE)
output_path = os.path.join(output_module_dir, OUTPUT_FILE)


# Check the UMI containing fastq file. Warning: Consumes a lot of memory; 1.5 times of the read file roughly
umis = dict()  # Load the input fastq to a dictionary
with open(read_with_umi, "r") as handle:  # Open the file
    for record in SeqIO.parse(handle, "fastq"):  # Read the fasta record by a function in Biopython package
        separated = record.id.split("____")
        if len(separated) == 2:  # In case there is no UMI for a given entry
            umis[separated[0]] = separated[1]  # The identifier is the key and UMI is the value


# Read the identifier, process and write to a new one
with open(input_path, "r") as handle:  # Read the fasta record by a function in Biopython package
    with open(output_path, "w") as output_handle:  # Open a file for output
        for record in SeqIO.parse(handle, "fasta"):
            if "____" not in record.id:  # Some entries have already UMI information
                record.id += "____" + umis[record.id]  # Add the info checking the dictionary
            output_handle.write(f">{record.id}\n{record.seq}\n")  # Write the new form of the entry
# Todo: Possible improvement: check if all entries preserved by looking the number of lines.


# End of the script
