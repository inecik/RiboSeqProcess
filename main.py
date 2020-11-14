"""
Analysis of Paired End Sequencing of Disome Profiling

The pipeline works only for a pair of reads.
If there is more than o sample or replicates, you have to run the script multiple times.
"""


import os
import subprocess
import sys
from shutil import which
from module_supplementary.common_functions import bcolors as c
import joblib


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


# Obtain the root dir for scripts to run properly
scripts_directory = os.path.dirname(__file__)  # Where this package is





# Preprocessing Module
# ____________________

print(f"{c.HEADER}Preprocessing Module.{c.ENDC}")
preprocess_inputs = {'read_1': sys.argv[1], 'read_2': sys.argv[2]}  # Command line inputs for the reads
# spr_1 = subprocess.run((
#     f"{which('python3')} "  # Define which python installation to use
#     f"{os.path.join(scripts_directory, 'module_preprocessing/cutadapt_umitools.py')} "  # Relative script dir
#     f"{preprocess_inputs['read_1']} "  # sys.argv[1]
#     f"{preprocess_inputs['read_2']} "  # sys.argv[2]
#     f"{OUTP_MAIN} "
#     f"{TEMP_MAIN}"
# ), shell=True)
# if spr_1.returncode != 0: sys.exit(f"{c.FAIL}Error in cutadapt_umitools.py: Exiting.{c.ENDC}")
print(f"{c.OKCYAN}Preprocessing Module was successfully completed.{c.ENDC}")





# Cleanup Module
# ______________

print(f"{c.HEADER}Clean-up Module.{c.ENDC}")
cleanup_inputs = joblib.load(os.path.join(TEMP_MAIN, ".module_preprocessing_paths.joblib"))

# Genome indexes
# spr_2 = subprocess.run((
#     f"{which('python3')} "  # Define which python installation to use
#     f"{os.path.join(scripts_directory, 'module_cleanup/database_rnaremove_bowtie2.py')} "  # Relative script dir
#     f"{TEMP_MAIN} "  # sys.argv[1]
# ), shell=True)
# if spr_2.returncode != 0: sys.exit(f"{c.FAIL}Error in database_rnaremove_bowtie2.py: Exiting.{c.ENDC}")
print(f"{c.OKCYAN}Genome indexing for rRNA removal was successfully completed.{c.ENDC}")

# Alignment
# spr_3 = subprocess.run((
#     f"{which('python3')} "  # Define which python installation to use
#     f"{os.path.join(scripts_directory, 'module_cleanup/bowtie2_rnaremove.py')} " 
#     f"{cleanup_inputs[0]} "  # sys.argv[1] Read 1
#     f"{cleanup_inputs[1]} "  # sys.argv[2] Read 2
#     f"{OUTP_MAIN} "  # sys.argv[3]
#     f"{TEMP_MAIN}"  # sys.argv[3]
# ), shell=True)
# if spr_3.returncode != 0: sys.exit(f"{c.FAIL}Error in bowtie2_rnaremove.py: Exiting.{c.ENDC}")
print(f"{c.OKCYAN}rRNA removal was successfully completed.{c.ENDC}")





# Link Pairing module
# ___________________

print(f"{c.HEADER}Link Pairing Module.{c.ENDC}")
linkpair_inputs = joblib.load(os.path.join(TEMP_MAIN, ".module_cleanup_paths.joblib"))

# Genome indexing
# spr_4 = subprocess.run((
#     f"{which('python3')} "  # Define which python installation to use
#     f"{os.path.join(scripts_directory, 'module_linkpairs/database_transcriptome_bowtie2.py')} "  # Relative script dir
#     f"{TEMP_MAIN} "  # sys.argv[1]
# ), shell=True)
# if spr_4.returncode != 0: sys.exit(f"{c.FAIL}Error in database_transcriptome_bowtie2.py: Exiting.{c.ENDC}")
print(f"{c.OKCYAN}Genome indexing for linking pairs was successfully completed.{c.ENDC}")

# Alignment
# spr_5 = subprocess.run((
#     f"{which('python3')} "  # Define which python installation to use
#     f"{os.path.join(scripts_directory, 'module_linkpairs/bowtie2_prealignment.py')} " 
#     f"{linkpair_inputs['read_1']} "  # sys.argv[1] 
#     f"{linkpair_inputs['read_2']} "  # sys.argv[2]
#     f"{OUTP_MAIN} "  # sys.argv[3]
#     f"{TEMP_MAIN}"  # sys.argv[3]
# ), shell=True)
# if spr_5.returncode != 0: sys.exit(f"{c.FAIL}Error in bowtie2_prealignment.py: Exiting.{c.ENDC}")
print(f"{c.OKCYAN}Alignment for linking pairs was successfully completed.{c.ENDC}")

# Create fasta
# spr_6 = subprocess.run((
#     f"{which('python3')} "  # Define which python installation to use
#     f"{os.path.join(scripts_directory, 'module_linkpairs/sam_processor.py')} " 
#     f"{OUTP_MAIN} "  # sys.argv[3]
#     f"{TEMP_MAIN}"  # sys.argv[3]
# ), shell=True)
# if spr_6.returncode != 0: sys.exit(f"{c.FAIL}Error in sam_processor.py: Exiting.{c.ENDC}")
print(f"{c.OKCYAN}Creation of fasta with linked pairs was successfully completed.{c.ENDC}")





# Genome alignment module
# _______________________

print(f"{c.HEADER}Genome Alignment Module.{c.ENDC}")
genomealignment_input = joblib.load(os.path.join(TEMP_MAIN, ".module_linkpair_paths.joblib"))

# Genome indexing
# spr_7 = subprocess.run((
#     f"{which('python3')} "  # Define which python installation to use
#     f"{os.path.join(scripts_directory, 'module_genomealignment/star_genome_index.py')} "  # Relative script dir
#     f"{TEMP_MAIN}"  # sys.argv[1]
# ), shell=True)
# if spr_7.returncode != 0: sys.exit(f"{c.FAIL}Error in star_genome_index.py: Exiting.{c.ENDC}")
print(f"{c.OKCYAN}Human genome indexing was successfully completed.{c.ENDC}")

# Alignment
spr_8 = subprocess.run((
    f"{which('python3')} "  # Define which python installation to use
    f"{os.path.join(scripts_directory, 'module_genomealignment/star_alignment_single.py')} " 
    f"{genomealignment_input} "  # sys.argv[1] 
    f"{OUTP_MAIN} "  # sys.argv[2]
    f"{TEMP_MAIN}"  # sys.argv[3]
), shell=True)
if spr_8.returncode != 0: sys.exit(f"{c.FAIL}Error in star_alignment_single.py: Exiting.{c.ENDC}")
print(f"{c.OKCYAN}Human genome alignment was successfully completed.{c.ENDC}")

# UMI deduplication
spr_9 = subprocess.run((
    f"{which('python3')} "  # Define which python installation to use
    f"{os.path.join(scripts_directory, 'module_genomealignment/umitools_dedup.py')} "  # Relative script dir
    f"{OUTP_MAIN}"
), shell=True)
if spr_9.returncode != 0: sys.exit(f"{c.FAIL}Error in umitools_dedup.py: Exiting.{c.ENDC}")
print(f"{c.OKCYAN}UMI deduplication was successfully completed.{c.ENDC}")


# Run assignment
