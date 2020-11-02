"""
This script is to find out the genes genomic positions.
"""


import sys
import os
import joblib


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
GFF3_FILE = "Homo_sapiens.GRCh38.96.chr.gff3"
OUTPUT_DATA_REPO = "id_mapping"
TEMP_DATA_REPO = "id_mapping"
OUTPUT_FILE = "genes_exones_genome_positions.joblib"
INPUT_FILE_1 = "transcript_geneelements.joblib"  # From genome_elements_mapping.py
INPUT_FILE_2 = "gene_transcript.joblib"  # From genome_elements_mapping.py


# Load the mapping data
data_repo_dir = os.path.join(output_dir, OUTPUT_DATA_REPO)
te_map = joblib.load(os.path.join(output_dir, data_repo_dir, INPUT_FILE_1))
gt_map = joblib.load(os.path.join(output_dir, data_repo_dir, INPUT_FILE_2))


gene_pos_dict = dict()  # Exons only
for gene in gt_map.keys():
    exons = [j for i in gt_map[gene] for j in te_map[i] if j[1] == 'exon']  # Loses 5' and 3' UTRs
    assert len(set([e[4] for e in exons])) == 1  # Check if all have the same orientation
    assert len(set([e[0] for e in exons])) == 1  # Check if all have the same chromosome
    positions = list()
    for e in exons:
        positions.extend(list(range(int(e[2]), int(e[3]))))  # Get the exon ranges
    positions = sorted(list(set(positions)))  # Drop overlapping positions, sort
    gene_pos_dict[gene] = (exons[0][0], positions)  # Chromosome and positions

# Write the results
print("Results is being written as a joblib object.")
output_path = os.path.join(output_dir, data_repo_dir, OUTPUT_FILE)
joblib.dump(gene_pos_dict, output_path)
print(f"Done: {output_path}")


# End of the script
