"""
This script creates two fastq file from a fasta file.
One with best possible sequencing quality and the other is worst possible sequencing quality.
It is to test whether STAR cares the sequencing quality when alignment.
"""




import os
import re
import sys
from Bio import SeqIO
import subprocess
from multiprocessing import cpu_count
from shutil import which


# Constants
best_chr = 'I'
worst_chr = '!'


# Operations for working environment and file name related operations
fasta_file = sys.argv[1]  # Important: Only file with extension ".fasta"
fasta_best = re.search(r"(.*)\.fasta$", fasta_file).group(1) + "_best.fastq"
fasta_worst = re.search(r"(.*)\.fasta$", fasta_file).group(1) + "_worst.fastq"


with open(fasta_file, "r") as handle_read:  # Open the file
    with open(fasta_best, "w") as handle_best, open(fasta_worst, "w") as handle_worst:
        for record in SeqIO.parse(handle_read, "fasta"):  # Read the fasta record by Biopython function
            handle_best.write("\n".join(["@" + record.id, str(record.seq),
                                         '+', best_chr * len(record.seq)]) + "\n")
            handle_worst.write("\n".join(["@" + record.id, str(record.seq),
                                          '+', worst_chr * len(record.seq)]) + "\n")


# End of the script
