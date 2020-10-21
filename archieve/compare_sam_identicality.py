"""
This script is to compare the results of two bam files on STAR.
It is a test to understand whether two SAM files are identical.
"""


import os
import re
import sys
import pysam
from Bio import SeqIO


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# Constants
sam1 = sys.argv[1]
sam2 = sys.argv[2]


aligned = dict()
counter = 0
with pysam.AlignmentFile(sam1, "r") as handle1:  # Open sam files
    sam_iterator = handle1.fetch()  # Get the iterator
    for e in sam_iterator:  # Iterate over the entries
        counter += 1
        if e.query_name in aligned:
            raise Exception("Multiple entry for the same identifier")
        else:
            aligned[e.query_name] = (e.reference_name, e.reference_start)  # Save the entry
print(f"Number of entries in first sam: {counter}")

counter = 0
with pysam.AlignmentFile(sam2, "r") as handle2:  # Open sam files
    sam_iterator = handle2.fetch()  # Get the iterator
    for e in sam_iterator:  # Iterate over the entries
        counter += 1
        try:
            assert aligned[e.query_name] == (e.reference_name, e.reference_start)
        except:
            print(f"The entry is different: {e.query_name}")
            print(f"Alignment in first sam: {aligned[e.query_name][0]}, {aligned[e.query_name][1]}")
            print(f"Alignment in second sam: {e.reference_name}, {e.reference_start}")
print(f"Number of entries in second sam: {counter}")


# End of the script
