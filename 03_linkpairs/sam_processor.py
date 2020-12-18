"""
The script parses SAM file to output pair's target and range.
bowtie2_prealignment.py should be run before this script.
It looks bowtie2_prealignment folder in the directory where the script is called.
The output folder will be mentioned bowtie2_prealignment folder.
"""


import os
import re
import sys
import joblib

import pysam
from Bio import SeqIO


# CONSTANTS
DATA_REPO = "03_linkpairs"  # Name of the database containing folder
SAM_FILE = "prealignment.sam"  # Check bowtie2_prealignment.py
OUTPUT_NAME = "footprints.fasta"
TRANSC_FASTA = "03_linkpairs/Homo_sapiens.GRCh38.cdna.all.fa"  # Check database_transc._bowtie2.py


# Operations for working environment and file name related operations
output_dir = sys.argv[1]
temp_dir = sys.argv[2]
data_repo_dir = os.path.join(output_dir, DATA_REPO)
sam_path = os.path.join(data_repo_dir, SAM_FILE)
fasta_transcriptome = os.path.join(temp_dir, TRANSC_FASTA)
output_fasta = os.path.join(data_repo_dir, OUTPUT_NAME)


# Take the transcript sequences into random access memory
print("Take the transcript sequences into random access memory.")
transcript_seqs = dict()  # Create a dictionary with transcript id as keys
with open(fasta_transcriptome, "r") as handle:  # Use previously filtered fasta file
    for record in SeqIO.parse(handle, "fasta"):
        transcript_seqs[record.id] = str(record.seq)  # Sequence as the values


# Processing the SAM file
print("Processing the SAM file.")
with open(output_fasta, "w") as output_handle:  # Open output fasta file
    popup_dict = dict()
    with pysam.AlignmentFile(sam_path, "r") as sam_handle:  # Open sam file
        sam_iterator = sam_handle.fetch()  # Get the iterator
        for e in sam_iterator:  # Iterate over the entries

            # Check if the entry is mapped, paired, and pair is mapped
            if (not e.is_unmapped and e.is_paired and e.is_proper_pair and not e.mate_is_unmapped):

                # Determine start and end positions, taking into account soft-clipping from the ends.
                qn = re.search(r"^[^_]*", e.query_name).group()  # Remove UMI info if exists
                is_entry_complete = False
                if not e.is_reverse and qn not in popup_dict:
                    popup_dict[qn] = [e.reference_start, None]
                elif not e.is_reverse:  # and qn in popup_dict
                    popup_dict[qn][0] = e.reference_start
                    is_entry_complete = True
                elif e.is_reverse and qn not in popup_dict:
                    popup_dict[qn] = [None, e.reference_end]
                elif e.is_reverse:  # and qn in popup_dict
                    popup_dict[qn][1] = e.reference_end
                    is_entry_complete = True
                else:
                    raise Exception("Unexpected entry!")

                # Print the result if the dict is fully completed.
                if is_entry_complete:
                    start_seq, end_seq = popup_dict.pop(qn)
                    # Different reference_names for pairs is impossible
                    fp = transcript_seqs[e.reference_name][start_seq: end_seq]
                    # Write down the identifier and sequence to fasta file
                    output_handle.write(f">{e.query_name}\n{fp}\n")


# Output the final file path to use in a pipeline
joblib.dump(output_fasta, os.path.join(output_dir, ".03_linkpairs.joblib"))


# End of the script
