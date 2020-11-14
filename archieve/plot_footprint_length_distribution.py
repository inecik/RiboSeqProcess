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
FIGURE_NAME = "footprint_length_distribution"


# Operations for working environment
input_path = sys.argv[1]
output_dir = sys.argv[2]
output_path_log = os.path.join(output_dir, FIGURE_NAME + "_log.pdf")
output_path_lin = os.path.join(output_dir, FIGURE_NAME + "_linear.pdf")


# Read fasta file and count the length of sequence
lengths = dict()  #
with open(input_path, "r") as handle:  # Read the fasta record by a function in Biopython package
    for record in SeqIO.parse(handle, "fasta"):
        l = len(record.seq)  # Get the length of the sequence
        if l in lengths:
            lengths[l] += 1  # Update to the dictionary if exists
        else:
            lengths[l] = 1  # Add to the dictionary if not exists


# Linear scale
seaborn.distplot(list(lengths.keys()), hist_kws={"weights":list(lengths.values())},
                 kde=False, norm_hist=True, bins=len(lengths), color="tomato")
plt.xlabel("Footprint Length")
plt.ylabel("Density")
seaborn.despine()  # Removes top and right frame
plt.savefig(output_path_log)  # Save to directory


# Log scale
seaborn.distplot(list(lengths.keys()), hist_kws={"weights":list(lengths.values())},
                 kde=False, norm_hist=True, bins=len(lengths), color="tomato")
plt.yscale('log')  # To see the distribution better
plt.xlabel("Footprint Length")
plt.ylabel("Density")
seaborn.despine()  # Removes top and right frame
plt.savefig(output_path_lin)  # Save to directory

# End of the script
