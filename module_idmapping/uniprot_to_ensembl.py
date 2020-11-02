"""
This script is to create a mapper between Uniprot and Ensembl.
"""


import sys
import os
import sys
import subprocess
import joblib
import pandas as pd


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# Inputs
output_dir = sys.argv[1]
temp_dir = sys.argv[2]


# CONSTANTS
MAPPING_DATA = "HUMAN_9606_idmapping.dat.gz"
OUTPUT_DATA_REPO = "id_mapping"
TEMP_DATA_REPO = "id_mapping"
OUTPUT_FILE = "uniprot_mapping.joblib"


# Create dir if not exist
data_repo_dir = os.path.join(output_dir, OUTPUT_DATA_REPO)
if not os.access(data_repo_dir, os.W_OK) or not os.path.isdir(data_repo_dir):  # Create directory if not exist
    os.mkdir(data_repo_dir)
# Create dir if not exist
temp_repo_dir = os.path.join(temp_dir, TEMP_DATA_REPO)
if not os.access(temp_repo_dir, os.W_OK) or not os.path.isdir(temp_repo_dir):  # Create directory if not exist
    os.mkdir(temp_repo_dir)


data_path_compressed = os.path.join(temp_dir, TEMP_DATA_REPO, MAPPING_DATA)
data_path = os.path.splitext(data_path_compressed)[0]
if not os.access(data_path, os.R_OK) or not os.path.isfile(data_path):
    subprocess.run((f"cd {temp_repo_dir}; curl -L -R -O ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/"
                    f"knowledgebase/idmapping/by_organism/{MAPPING_DATA}"), shell=True)
    subprocess.run(f"cd {temp_repo_dir}; gzip -d {data_path_compressed}", shell=True)


# Read the data file
with open(data_path, "r") as input_handle:
    data_raw = input_handle.readlines()
df = list()
for entry in data_raw:
    entry = entry.strip().split("\t")
    if entry[1] in ["Ensembl", "Ensembl_TRS", "Ensembl_PRO", "Gene_Name", "Gene_Synonym"]:
        df.append(entry)

# Write the results
print("Results is being written as a joblib object.")
output_path = os.path.join(output_dir, data_repo_dir, OUTPUT_FILE)
joblib.dump(df, output_path)
print(f"Done: {output_path}")


# End of the script
