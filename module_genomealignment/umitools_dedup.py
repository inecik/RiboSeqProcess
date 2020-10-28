"""
Star alignment for single end sequencing.

"""


import os
import re
import sys
import subprocess
from multiprocessing import cpu_count
from datetime import datetime
from shutil import which


# CONSTANTS
OUTPUT_STAR = "star_single_module"  # from star_alignment_single.py
PREFIX = "single_"  # from star_alignment_single.py
OUTPUT_PREFIX = "umi-deduplicated"


# Inputs
output_dir = sys.argv[1]  # For the directory for output files, as described in main.py
output_dir_module = os.path.join(output_dir, OUTPUT_STAR)
raw_sam = os.path.join(output_dir_module, PREFIX + "_Aligned.sortedByCoord.out.sam")
os.chdir(output_dir_module)   # Since everything will be output there


# Preprocessing
after_convert = os.path.splitext(raw_sam)[0] + '.bam'
subprocess.run((  # Convert sam file into bam file
    f"{which('samtools')} view -bS "  # Define which samtools installation to use
    f"{raw_sam} > "  # Input
    f"{after_convert}"  #
), shell=True)

after_sort = os.path.splitext(after_convert)[0] + '_sorted.bam'
subprocess.run((  # Sort bam file
    f"{which('samtools')} sort "  # Define which samtools installation to use
    f"{after_convert} "  # Input
    f"-o {after_sort}"  #
), shell=True)

subprocess.run((  # Index
    f"{which('samtools')} index "  # Define which samtools installation to use
    f"{after_sort}"  # Input
), shell=True)


# Deduplication of UMI
subprocess.run((  # Run deduplication
    f"{which('umi_tools')} dedup "  # Define which umi_tools installation to use
    f"-I {after_sort} "  # Input
    f"--output-stats={OUTPUT_PREFIX} "  
    f"-S {OUTPUT_PREFIX}.bam "
    f"> {OUTPUT_PREFIX}.log"
), shell=True)


# Results
subprocess.run((  # Convert to sam file
    f"{which('samtools')} view -h "  # Define which samtools installation to use
    f"-o {OUTPUT_PREFIX}.sam "  # Output
    f"{OUTPUT_PREFIX}.bam"  # Input
), shell=True)

subprocess.run((  # Convert to fasta file
    f"{which('samtools')} bam2fq "  # Define which samtools installation to use
    f"{OUTPUT_PREFIX}.bam "  # Input
    f"| "  # Pipe with the following
    f"{which('seqtk')} seq -A - > "
    f"{OUTPUT_PREFIX}.fasta"
), shell=True)


# End of the script
