"""
The script parses SAM file to output pair's target and range.
bowtie2_prealignment.py should be run before this script.
It looks bowtie2_prealignment folder in the directory where the script is called.
The output folder will be mentioned bowtie2_prealignment folder.
"""


import os
import pysam
from Bio import SeqIO


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# CONSTANTS
DATA_REPO = "bowtie2_prealignment"  # Name of the database containing folder
SAM_FILE = "bowtie2_prealignment.sam"  # Check bowtie2_prealignment.py
OUTPUT_NAME = "footprints.fasta"
TRANSC_FASTA = "bowtie2_transcriptome/Homo_sapiens.GRCh38.cdna.all_filtered.fa"  # Check database_transc._bowtie2.py


# Operations for working environment and file name related operations
running_directory = os.getcwd()
data_repo_dir = os.path.join(running_directory, DATA_REPO)
sam_path = os.path.join(data_repo_dir, SAM_FILE)
fasta_transcriptome = os.path.join(running_directory, TRANSC_FASTA)
output_fasta = os.path.join(data_repo_dir, OUTPUT_NAME)


# Take the transcript sequences into random access memory
transcript_seqs = dict()  # Create a dictionary with transcript id as keys
with open(fasta_transcriptome, "r") as handle:  # Use previously filtered fasta file
    for record in SeqIO.parse(handle, "fasta"):
        transcript_seqs[record.id] = str(record.seq)  # Sequence as the values


# Processing the SAM file
with pysam.AlignmentFile(sam_path, "r") as sam_handle:  # Open sam file
    with open(output_fasta, "w") as output_handle:  # Open output fasta file
        sam_iterator = sam_handle.fetch()  # Get the iterator
        for e in sam_iterator:  # Iterate over the entries
            if (not e.is_unmapped  # Check if the entry is mapped, paired, and pair is mapped
                    and e.is_paired
                    and e.is_proper_pair
                    and not e.mate_is_unmapped
                    and e.template_length > 0  # Take only the read at forward position. It is not necessarily read 1.
                    ):  # If the read satisfies all above
                # Get the entry reference name and fetch the associated sequence from transcriptome dictionary
                # Substring the reference by ranges [from Start to Start + Inferred length]
                fp = transcript_seqs[e.reference_name][e.reference_start: e.reference_start + e.template_length]

                # todo: copy and paste 5' and 3' end instead of taking directly from the reference genome
                # e.template_length > 0 ve e.template_length < 0 ile iki çifti al
                # UMI'lerden arındırarak bunları mate ilan et,

                # reverse transcriptase'ın verdiği rastgele şeyleri alıyor. şu andaki footprinter

                output_handle.write(f">{e.query_name}\n{fp}\n")  # Write down the identifier and sequence to fasta file


# End of the script
