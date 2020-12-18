"""
Star alignment for single end sequencing.

"""


import os
import sys
import subprocess
from shutil import which
import joblib

repository_name = "RiboSeqProcess"
sys.path.append(os.path.abspath(__file__).split(repository_name)[0] + repository_name)
from archieve.common_functions import *


check_exist_package("/usr/bin/samtools")  # use system samtools
check_exist_package("umi_tools")
check_exist_package("seqtk")


# CONSTANTS
OUTPUT_STAR = "04_genomealignment"  # from star_alignment_single.py
PREFIX = "paired_"  # from star_alignment_single.py
OUTPUT_PREFIX = "umi-deduplicated"


# Inputs
output_dir = sys.argv[1]  # For the directory for output files, as described in main.py
output_dir_module = os.path.join(output_dir, OUTPUT_STAR)
raw_bam = os.path.join(output_dir_module, PREFIX + "Aligned.sortedByCoord.out.bam")
os.chdir(output_dir_module)   # Since everything will be output there


# Preprocessing
print("Bam file is being sorted.")
after_sort = os.path.splitext(raw_bam)[0] + '_sorted.bam'
subprocess.run((  # Sort bam file
    f"{which('/usr/bin/samtools')} sort "  # Define which samtools installation to use
    f"{raw_bam} "  # Input
    f"-o {after_sort}"  #
), shell=True)

print("Bam file is being indexed.")
subprocess.run((  # Index
    f"{which('/usr/bin/samtools')} index "  # Define which samtools installation to use
    f"{after_sort}"  # Input
), shell=True)


# Deduplication of UMI
print("UMI deduplication process is started.")
subprocess.run((  # Run deduplication
    f"{which('umi_tools')} dedup "  # Define which umi_tools installation to use
    f"-I {after_sort} "  # Input
    f"--output-stats={OUTPUT_PREFIX} "
    f"--paired "
    f"--unpaired-reads discard "  
    f"-S {OUTPUT_PREFIX}.bam "
    f"> {OUTPUT_PREFIX}.log"
), shell=True)


# Results
print("Sam file is being created.")
subprocess.run((  # Convert to sam file
    f"{which('/usr/bin/samtools')} view -h "  # Define which samtools installation to use
    f"-o {OUTPUT_PREFIX}.sam "  # Output
    f"{OUTPUT_PREFIX}.bam"  # Input
), shell=True)


final_paths = {'sam': os.path.join(output_dir_module, f"{OUTPUT_PREFIX}.sam"),
               'bam': os.path.join(output_dir_module, f"{OUTPUT_PREFIX}.bam")}
joblib.dump(final_paths, os.path.join(output_dir, ".04_genomealignment.joblib"))


# End of the script
