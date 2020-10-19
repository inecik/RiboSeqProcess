"""
It is to look at the distribution of footprint length.
You should run this script where you have main.
It assumes there is a file called bowtie2_prealignment/footprints_UMI.fasta in running dir.
"""


import os
from Bio import SeqIO
import sys
import matplotlib.pyplot as plt
import seaborn

# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# CONSTANTS
FIGURE_MAIN = "figures"  # name of the figure containing folder
DATABASE_FASTA = "bowtie2_prealignment/footprints_UMI.fasta"  # Name of footprint fasta file
FIGURE_NAME = "footprint_length_distribution.pdf"


# Operations for working environment
running_directory = os.getcwd()
data_repo_dir = os.path.join(running_directory, FIGURE_MAIN)
if not os.access(data_repo_dir, os.W_OK) or not os.path.isdir(data_repo_dir):  # Create directory if not exist
    print("Data directory created")
    os.mkdir(data_repo_dir)
input_path = os.path.join(running_directory, DATABASE_FASTA)
output_path = os.path.join(data_repo_dir, FIGURE_NAME)


# Read fasta file and count the length of sequence
lengths = dict()  #
with open(input_path, "r") as handle:  # Read the fasta record by a function in Biopython package
    for record in SeqIO.parse(handle, "fasta"):
        l = len(record.seq)  # Get the length of the sequence
        if l in lengths:
            lengths[l] += 1  # Update to the dictionary if exists
        else:
            lengths[l] = 1  # Add to the dictionary if not exists


seaborn.distplot(list(lengths.keys()), hist_kws={"weights":list(lengths.values())},
                 kde=False, norm_hist=True, bins=len(lengths))
plt.yscale('log')  # To see the distribution better
plt.xlabel("Footprint Length")
plt.ylabel("Density")
seaborn.despine()  # Removes top and right frame
plt.savefig(output_path)  # Save to directory


# End of the script
