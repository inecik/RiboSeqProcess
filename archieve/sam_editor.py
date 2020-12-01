"""
Sam editor to test sam-tools with Ilia's version of writing UMI
"""

import os
import pysam



# Operations for working environment and file name related operations
read_sam = "/Users/kemalinecik/Desktop/single.sam"  # sys.argv[1]
write_sam = os.path.splitext(read_sam)[0] + "_edited.sam"

with pysam.AlignmentFile(read_sam, "r") as input_handle:  # Open output fasta file
    with pysam.AlignmentFile(write_sam, "w", template=input_handle) as output_handle:  # Open output fasta file
        sam_iterator = input_handle.fetch()  # Get the iterator
        for e in sam_iterator:  # Iterate over the entries
            e.query_name = "_".join(e.query_name.split("____"))
            output_handle.write(e)


# End of the script
