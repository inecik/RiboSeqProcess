"""
The script parses SAM file to output pair's target and range.
bowtie2_prealignment.py should be run before this script.
It looks bowtie2_prealignment folder in the directory where the script is called.
The output folder will be mentioned bowtie2_prealignment folder.
"""


import os
import re
import subprocess
from multiprocessing import cpu_count
from shutil import which
import pysam
from Bio import SeqIO


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# Operations for working environment
running_directory = os.getcwd()
data_repository = "bowtie2_prealignment"  # Name of the database containing folder
data_repo_dir = os.path.join(running_directory, data_repository)

# Todo: here sloppy
sam_path = os.path.join(data_repo_dir, "bowtie2_prealignment.sam")
fasta_transcriptome =  "/Users/kemalinecik/Documents/Programming/Ribo-seq-Analysis/bowtie2_transcriptome/Homo_sapiens.GRCh38.cdna.all_filtered.fa"
output_fasta = os.path.join(data_repo_dir, "footprints.fasta")


# Take the transcript sequences into random access memory
transcript_seqs = dict()
with open(fasta_transcriptome, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        transcript_seqs[record.id] = str(record.seq)


# UMI dataframe
# todo

with pysam.AlignmentFile(sam_path, "r") as sam_handle:  # Open sam file
    with open(output_fasta, "w") as output_handle:
        sam_iterator = sam_handle.fetch()  # Get the iterator
        for e in sam_iterator:
            if (not e.is_unmapped
                    and e.is_paired
                    and e.is_proper_pair
                    and not e.mate_is_unmapped
                    and e.template_length > 0
                    ):
                fp = transcript_seqs[e.reference_name][e.reference_start : e.reference_start + e.template_length]
                output_handle.write(f">{e.query_name}\n{fp}\n")

