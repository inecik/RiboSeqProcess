"""
Analysis of Paired End Sequencing of Disome Profiling

The working directory is the directory where the code is run.
The pipeline works only for a pair of reads.
If there is more than o sample or replicates, you have to run the script multiple times.

The output file is created at the working directory.
"""


import os
import subprocess
import sys
from shutil import which
import multiprocessing


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# Inputs
OUTP_MAIN = sys.argv[3]
TEMP_MAIN = sys.argv[4]
FIGU_MAIN = sys.argv[5]


# Obtain the root dir for scripts to run properly
scripts_directory = os.path.dirname(__file__)  # Where this package is


# Preprocessing Module
#
preprocess_inputs = {'read_1': sys.argv[1], 'read_2': sys.argv[2]}  # Command line inputs for the reads
subprocess.run((
    f"{which('python3')} "  # Define which python installation to use
    f"{os.path.join(scripts_directory, 'module_preprocessing/cutadapt_umitools.py')} "  # Relative script dir
    f"{preprocess_inputs['read_1']} "  # sys.argv[1] 
    f"{preprocess_inputs['read_2']} "  # sys.argv[2]
    f"{OUTP_MAIN}"
), shell=True)


# Cleanup Module
#
cleanup_inputs = {'read_1': "read_1_no-adapt_umi-aware.fastq.gz", 'read_2': "read_2_cutadapt_umi_aware.fastq.gz"}

# Genome indexes
subprocess.run((
    f"{which('python3')} "  # Define which python installation to use
    f"{os.path.join(scripts_directory, 'module_cleanup/database_rnaremove_bowtie2.py')} "  # Relative script dir
    f"{TEMP_MAIN} "  # sys.argv[1]
), shell=True)

# Alignment
subprocess.run((
    f"{which('python3')} "  # Define which python installation to use
    f"{os.path.join(scripts_directory, 'module_cleanup/bowtie2_rnaremove.py')} " 
    f"{cleanup_inputs['read_1']} "  # sys.argv[1] 
    f"{cleanup_inputs['read_2']} "  # sys.argv[2]
    f"{OUTP_MAIN} "  # sys.argv[3]
    f"{TEMP_MAIN}"  # sys.argv[3]
), shell=True)


# Link Pairing module
#
linkpair_inputs = {'read_1': "read_1_norRNA.fastq", 'read_2': "read_2_norRNA.fastq"}

# Genome indexing
subprocess.run((
    f"{which('python3')} "  # Define which python installation to use
    f"{os.path.join(scripts_directory, 'module_linkpairs/database_transcriptome_bowtie2.py')} "  # Relative script dir
    f"{TEMP_MAIN} "  # sys.argv[1]
), shell=True)

# Alignment
subprocess.run((
    f"{which('python3')} "  # Define which python installation to use
    f"{os.path.join(scripts_directory, 'module_linkpairs/bowtie2_prealignment.py')} " 
    f"{linkpair_inputs['read_1']} "  # sys.argv[1] 
    f"{linkpair_inputs['read_2']} "  # sys.argv[2]
    f"{OUTP_MAIN} "  # sys.argv[3]
    f"{TEMP_MAIN}"  # sys.argv[3]
), shell=True)

# Re-create fasta
subprocess.run((
    f"{which('python3')} "  # Define which python installation to use
    f"{os.path.join(scripts_directory, 'module_linkpairs/sam_processor.py')} " 
    f"{OUTP_MAIN} "  # sys.argv[3]
    f"{TEMP_MAIN}"  # sys.argv[3]
), shell=True)


# Genome alignment module
#
genomealignment_inputs = {"read": "bowtie2_prealignment/footprints.fasta"}

# Genome indexing
subprocess.run((
    f"{which('python3')} "  # Define which python installation to use
    f"{os.path.join(scripts_directory, 'module_genomealignment/star_genome_index.py')} "  # Relative script dir
    f"{TEMP_MAIN}"  # sys.argv[1]
), shell=True)

# Alignment
subprocess.run((
    f"{which('python3')} "  # Define which python installation to use
    f"{os.path.join(scripts_directory, 'module_genomealignment/star_alignment_single.py')} " 
    f"{genomealignment_inputs['read']} "  # sys.argv[1] 
    f"{OUTP_MAIN} "  # sys.argv[2]
    f"{TEMP_MAIN}"  # sys.argv[3]
), shell=True)

# UMI deduplication
subprocess.run((
    f"{which('python3')} "  # Define which python installation to use
    f"{os.path.join(scripts_directory, 'module_genomealignment/umitools_dedup.py')} "  # Relative script dir
    f"{OUTP_MAIN}"
), shell=True)

# Run assignment
