"""

"""

import os
import re
import sys
from Bio import SeqIO
import gzip

inp = sys.argv[1]
length = int(sys.argv[2])
output = os.path.splitext(inp)[0] + f"_first_{length}.fastq.gz"

# Take the query names into the RAM for Read 1
print("Take the query names into the RAM for Read 1.")
read_1_set = list()
with gzip.open(inp, "rt") as handle:  # Open the gz file
    for ind, record in enumerate(SeqIO.parse(handle, "fastq")):  # Read the fasta record by a function in Biopython package
        read_1_set.append(record)
        if ind + 1 >= length:
            break

print("Write down the new file")
with gzip.open(output, 'wt') as handle_1:
    for record in read_1_set:
        SeqIO.write(record, handle_1, "fastq")



