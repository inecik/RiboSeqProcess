"""

"""


import sys
import os
import subprocess
from shutil import which

repository_name = "RiboSeqProcess"
sys.path.append(os.path.abspath(__file__).split(repository_name)[0] + repository_name)
from archieve.common_functions import *

# Inputs
sam_file = sys.argv[1]
output_dir = sys.argv[2]
temp_dir = sys.argv[3]


# CONSTANTS
GFF3_FILE = "Homo_sapiens.GRCh38.96.chr.gff3.gz"
OUTPUT_DATA_REPO = "05_julia_assignment"
TEMP_DATA_REPO = "05_julia_assignment"


# Create dir if not exist
data_repo_dir = create_dir(output_dir, OUTPUT_DATA_REPO)
temp_repo_dir = create_dir(temp_dir, TEMP_DATA_REPO)


# Create or load SQL database
gff3_uncomp = os.path.splitext(GFF3_FILE)[0]
gff_path = os.path.join(temp_dir, TEMP_DATA_REPO, gff3_uncomp)
corrected_gff = os.path.splitext(gff_path)[0] + "_renamed_duplicate_gene_names.gff3"
if not os.access(corrected_gff, os.R_OK) or not os.path.isfile(corrected_gff):
    if not os.access(gff_path, os.R_OK) or not os.path.isfile(gff_path):
        subprocess.run((f"cd {temp_repo_dir}; curl -L -R -O ftp://ftp.ensembl.org/pub/release-96/gff3/"
                       "homo_sapiens/Homo_sapiens.GRCh38.96.chr.gff3.gz"), shell=True)
        subprocess.run(f"cd {temp_repo_dir}; gzip -d Homo_sapiens.GRCh38.96.chr.gff3.gz", shell=True)
    subprocess.run((f"cd {os.path.dirname(__file__)}; "
                   f"{which('python3')} gff3_rename_duplicated_genes.py {gff_path}"), shell=True)


# Run Ilia's Julia script
julia_path = os.path.join(os.path.dirname(__file__), "02_assign_count_streaming_11.10.19_modified.jl")
subprocess.run((
    f"{which('julia')} {julia_path} "  # Which Julia installation to use and the script
    f"-g {corrected_gff} "  # Gff3 file. Removed of duplicated gene names
    "-a 3 "  # Assignment from 5'
    # "-u "  # Inherited from Mati. Removed because umi-tool deduplication is already done.
    f"-o {data_repo_dir} "  # Output file
    f"{sam_file}"
), shell=True)


# End of the script
