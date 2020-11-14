"""
Umitools for paired end sequencing, runs for both reads.

It is the same as the umi_tools extract in cutadapt.py but it also syncronizes the fastq files.
"""


import os
import sys
import gzip
from Bio import SeqIO
import re

sys.path.append(os.path.abspath(__file__).split("Ribo-seq-Analysis")[0] + "Ribo-seq-Analysis")

# SETTINGS
PATTERN_read1 = "^(?P<umi_1>.{2})(.*)$"  # Regex to get UMI from first 2 in read1
PATTERN_read2 = "^(.*)(?P<umi_2>.{5})(?P<discard_2>.{5})$"  # # Regex to get rid of barcode and take 5 nt UMI from read 2


# Inputs
output_dir = sys.argv[1]  # Command line input 3 for the directory for output files, as described in main.py
OUTPUT_DATA_REPO = "preprocessing_module"  # from cutadapt.py
output_dir_module = os.path.join(output_dir, OUTPUT_DATA_REPO)  # from cutadapt.py
os.chdir(output_dir_module)  # Since everything will be output there

NAMES = ["read_1", "read_2"]  # from cutadapt.py
temp_paths = [f"{NAMES[0]}_cutadapt_temp.fastq.gz", f"{NAMES[1]}_cutadapt_temp.fastq.gz"]  # from cutadapt.py
temp_paths = [os.path.join(output_dir_module, i) for i in temp_paths]  # from cutadapt.py

# UMI extracted outputs
final_paths = [f"{i}_no-adapt_umi-aware_sync.fastq.gz" for i in NAMES]
final_paths = [os.path.join(output_dir_module, i) for i in final_paths]


# Take the query names into the RAM for Read 1
print("Take the query names into the RAM for Read 1.")
read_1_set = set()
read_1_dict = dict()
with gzip.open(temp_paths[0], "rt") as handle:  # Open the gz file
    for record in SeqIO.parse(handle, "fastq"):  # Read the fasta record by a function in Biopython package
        read_1_set.add(record.name)
        read_1_dict[record.name] = record

print("Take the query names into the RAM for Read 2.")
read_2_set = set()
read_2_dict = dict()
with gzip.open(temp_paths[1], "rt") as handle:  # Open the gz file
    for record in SeqIO.parse(handle, "fastq"):  # Read the fasta record by a function in Biopython package
        read_2_set.add(record.name)
        read_2_dict[record.name] = record

print("Find the intersections")
intersect_records = read_1_set.intersection(read_2_set)

print("Extract UMIs & Sync fastq files")
with gzip.open(final_paths[0], 'wt') as handle_1, gzip.open(final_paths[1], 'wt') as handle_2:
    for qn in intersect_records:
        # For read 1
        read_1 = read_1_dict[qn]
        umi_1, trimmed_read_1 = re.search(PATTERN_read1, read_1.seq._data).groups()
        read_1 = read_1[2:]

        # For read 2
        read_2 = read_2_dict[qn]
        trimmed_read_2, umi_2, barcode_2 = re.search(PATTERN_read2, read_2.seq._data).groups()
        read_2 = read_2[:-10]

        # UMI infos to the title
        read_1.name = read_1.id = read_1.description = qn + "_" + umi_1 + umi_2
        read_2.name = read_2.id = read_2.description = qn + "_" + umi_1 + umi_2

        # Write the result
        SeqIO.write(read_1, handle_1, "fastq")
        SeqIO.write(read_2, handle_2, "fastq")


# End of the script
