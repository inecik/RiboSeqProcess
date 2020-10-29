"""
This script reads GERP scores and single-end SAM file.
Finds out the conservation score of the reads and genomic positions.
It requires 150 GBs of RAM
"""


import os
import sys
import pyBigWig
import joblib
import pandas as pd

sys.path.append(os.path.abspath(__file__).split("Ribo-seq-Analysis")[0] + "Ribo-seq-Analysis")
from archieve.common_functions import progressBarForTerminal

# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# Inputs
output_dir = "/home/kai/KEMALINECIK/out/OUTPUT"
temp_dir = "/home/kai/KEMALINECIK/out/TEMP"


# CONSTANTS
GERP_BIGWIG = "gerp_conservation_scores.homo_sapiens.GRCh38.bw"
OUTPUT_DATA_REPO = "gerp_genes"
TEMP_DATA_REPO = "gerp"
GENE_POSITIONS = "genes_genome_positions.joblib"
OUTPUT_FILE = "genes_positions_cons-scores.joblib"


# Create dir if not exist
data_repo_dir = os.path.join(output_dir, OUTPUT_DATA_REPO)
if not os.access(data_repo_dir, os.W_OK) or not os.path.isdir(data_repo_dir):  # Create directory if not exist
    os.mkdir(data_repo_dir)
# Create dir if not exist
temp_repo_dir = os.path.join(temp_dir, TEMP_DATA_REPO)
if not os.access(temp_repo_dir, os.W_OK) or not os.path.isdir(temp_repo_dir):  # Create directory if not exist
    os.mkdir(temp_repo_dir)


# Read genome positions
gene_pos_dict_path = os.path.join(output_dir, data_repo_dir, GENE_POSITIONS)
gene_pos_dict = joblib.load(gene_pos_dict_path)
gerp_bigwig_path = os.path.join(temp_repo_dir, GERP_BIGWIG)


# Take all conservation file into memory
cons = dict()
with pyBigWig.open(gerp_bigwig_path, 'r') as bw:
    print("Reading conservation scores of chromosome:")
    for chr in bw.chroms():
        print(chr)
        cons[chr] = bw.values(chr, 0, bw.chroms()[chr])
print("Score reading completed")


# Fetch the score and crete a dataframe to store the data
gerp_bigwig_path = os.path.join(temp_repo_dir, GERP_BIGWIG)
df_genes_out = pd.DataFrame(columns=["query_name", "chromosome", "positions", "scores"])
reporter = 0
with pyBigWig.open(gerp_bigwig_path, "r") as bw:
    for gene in gene_pos_dict.keys():  # Iterate over the genes
        reporter += 1
        gene_pos = gene_pos_dict[gene]
        scores = [cons[gene_pos[0]][i] for i in gene_pos[1]]
        df_genes_out.append({
            "query_name": gene,
            "chromosome": gene_pos[0],
            "positions": gene_pos[1],
            "scores": scores,
        }, ignore_index=True)
        if reporter == gene_pos_dict-1 or reporter % 10000 == 0:
            progressBarForTerminal(reporter, gene_pos_dict-1)


# Write the results
output_path = os.path.join(data_repo_dir, OUTPUT_FILE)
joblib.dump(df_genes_out, output_path)


# End of the script
