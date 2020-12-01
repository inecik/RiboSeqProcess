"""
The script downloads and filters non-coding RNA data for Bowtie2 alignment tool.
It filters out everything other than rRNA, tRNA
Then, it runs indexing tool of Bowtie2.
The output folder will be created to the directory where the script is called.
"""


import os
import subprocess
import joblib
import sys
from multiprocessing import cpu_count

repository_name = "RiboSeqProcess"
sys.path.append(os.path.abspath(__file__).split(repository_name)[0] + repository_name)
from archieve.common_functions import *


# Check if necessary packages were installed.
check_exist_package("bowtie2-build")


# CONSTANTS
MODULE_TEMP = "02_cleanup"  # name of the database containing folder
INDEX_BASE = "homo_sapiens_rrna_trna"  # The basename of the index files to write


# Operations for working environment
temp_dir = sys.argv[1]  # Check main.py for it
temp_module_dir = create_dir(temp_dir, MODULE_TEMP)  # Create directory if not exist


# Search https://rnacentral.org to the following
# TAXONOMY:"9606" AND (so_rna_type_name:"RRNA" OR so_rna_type_name:"TRNA" ) AND (expert_db:"Ensembl" OR expert_db:"SILVA" OR expert_db:"NONCODE" OR expert_db:"Rfam" OR expert_db:"PDBe" OR expert_db:"RefSeq" OR expert_db:"GeneCards" ) AND (rna_type:"rRNA" OR rna_type:"tRNA")


# Path of the manual curation of human rRNAs
non_coding_fasta_name = "human_rRNAs_tRNAs.fasta"
non_coding_fasta = os.path.join(os.path.abspath(os.path.dirname(__file__)), non_coding_fasta_name)
check_exist_file(non_coding_fasta)


# Run Bowtie2_build function in Bowtie2 module.
bowtie2_index = create_dir(temp_module_dir, "index_ncrna")
metadata_index_path = os.path.join(temp_module_dir, ".02_cleanup_metadata.joblib")

try:
    # Check if the index file
    assert os.path.isfile(metadata_index_path) and os.access(metadata_index_path, os.R_OK)
    metadata_index_previously = joblib.load(metadata_index_path)
    metadata_index = get_files_metadata(bowtie2_index)
    assert metadata_index_previously == metadata_index
    print("Index files for rRNA removal are already present.")

except:
    print("Index files for rRNA removal are now being created.")
    subprocess.run((
        f"cd {bowtie2_index};"  # Change the directory to the index directory
        f"{which('bowtie2-build')} "  # Name of the function
        f"--threads {cpu_count()} "  # Number of threads to be used.
        f"{non_coding_fasta} "  # Input file. -f is to indicate the file is in fasta format
        f"{INDEX_BASE} "  # The basename of the index files to write
        "> report_indexing.log"
    ), shell=True)

    # Write the file info
    metadata_index = get_files_metadata(bowtie2_index)
    joblib.dump(metadata_index, metadata_index_path)


# End of the script
