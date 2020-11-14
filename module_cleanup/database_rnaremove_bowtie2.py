"""
The script downloads and filters non-coding RNA data for Bowtie2 alignment tool.
It filters out everything other than rRNA, tRNA
Then, it runs indexing tool of Bowtie2.
The output folder will be created to the directory where the script is called.
"""


import os
import subprocess
from shutil import which
import sys
from multiprocessing import cpu_count

sys.path.append(os.path.abspath(__file__).split("Ribo-seq-Analysis")[0] + "Ribo-seq-Analysis")
from module_supplementary.common_functions import bcolors as c


# Check if necessary packages were installed.
if not which('bowtie2-build'):
    sys.exit(f"{c.FAIL}Bowtie2 package should be installed.")
else:
    print("Bowtie2 installation is found.")


# CONSTANTS
MODULE_TEMP = "bowtie2_rrna-trna"  # name of the database containing folder
INDEX_BASE = "homo_sapiens_rrna_trna"  # The basename of the index files to write


# Operations for working environment
temp_dir = sys.argv[1]  # Check main.py for it
temp_module_dir = os.path.join(temp_dir, MODULE_TEMP)
if not os.access(temp_module_dir, os.W_OK) or not os.path.isdir(temp_module_dir):  # Create directory if not exist
    os.mkdir(temp_module_dir)
    print("Temp directory created.")


# Search https://rnacentral.org to the following
# TAXONOMY:"9606" AND (so_rna_type_name:"RRNA" OR so_rna_type_name:"TRNA" ) AND (expert_db:"Ensembl" OR expert_db:"SILVA" OR expert_db:"NONCODE" OR expert_db:"Rfam" OR expert_db:"PDBe" OR expert_db:"RefSeq" OR expert_db:"GeneCards" )


# Path of the manual curation of human rRNAs
non_coding_fasta_name = "human_rRNAs_tRNAs.fasta"
non_coding_fasta = os.path.join(os.path.abspath(os.path.dirname(__file__)), non_coding_fasta_name)
if not os.access(temp_module_dir, os.R_OK) or not os.path.isfile(non_coding_fasta):
    sys.exit(f"{c.FAIL}Manually downloaded ncRNA fasta file is not found.")


# Run Bowtie2_build function in Bowtie2 module.
print("Genome indexes for rRNA removal is now being created.")
bowtie2_index = os.path.join(temp_module_dir, "bowtie2_index")
if not os.access(bowtie2_index, os.W_OK) or not os.path.isdir(bowtie2_index):  # Create directory if not exist
    print("Genome indexing directory is created under temp dir.")
    os.mkdir(bowtie2_index)
subprocess.run((
    f"cd {bowtie2_index};"  # Change the directory to the index directory
    f"{which('bowtie2-build')} "  # Name of the function
    f"--threads {cpu_count()} "  # Number of threads to be used.
    f"{non_coding_fasta} "  # Input file. -f is to indicate the file is in fasta format
    f"{INDEX_BASE} "  # The basename of the index files to write
    "> report_indexing.log"
), shell=True)


# End of the script
