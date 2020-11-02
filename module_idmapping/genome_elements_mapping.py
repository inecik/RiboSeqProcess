"""
This script is to have mapping tool transcripts to transcripts elements and genes to transcripts
"""


import sys
import os
import subprocess
import joblib
from itertools import compress


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
OUTPUT_FILE_1 = "transcript_geneelements.joblib"
OUTPUT_FILE_2 = "gene_transcript.joblib"


# Create dir if not exist
data_repo_dir = os.path.join(output_dir, OUTPUT_DATA_REPO)
if not os.access(data_repo_dir, os.W_OK) or not os.path.isdir(data_repo_dir):  # Create directory if not exist
    os.mkdir(data_repo_dir)
# Create dir if not exist
temp_repo_dir = os.path.join(temp_dir, TEMP_DATA_REPO)
if not os.access(temp_repo_dir, os.W_OK) or not os.path.isdir(temp_repo_dir):  # Create directory if not exist
    os.mkdir(temp_repo_dir)


# Check the presence of the database
gff_path = os.path.join(temp_dir, TEMP_DATA_REPO, GFF3_FILE)
if not os.access(gff_path, os.R_OK) or not os.path.isfile(gff_path):
    subprocess.run((f"cd {temp_repo_dir}; curl -L -R -O ftp://ftp.ensembl.org/pub/release-96/gff3/"
                   "homo_sapiens/Homo_sapiens.GRCh38.96.chr.gff3.gz"), shell=True)
    subprocess.run(f"cd {temp_repo_dir}; gzip -d Homo_sapiens.GRCh38.96.chr.gff3.gz", shell=True)


with open(gff_path, "r") as gff_handle:
    gff_raw = gff_handle.readlines()
te_map = dict()  # Transcript genomic element map
gt_map = dict()  # Gene transcript map
for entry in gff_raw:
    if not entry.startswith('#'):
        entry = entry.strip().split('\t')
        attributes = dict([i.split('=') for i in entry[8].split(';')])
        if entry[2] in ["mRNA", "exon", "five_prime_UTR", "three_prime_UTR", "CDS"]:
            parent = attributes['Parent']
            if parent.startswith('transcript:'):
                transcript_id = parent.split(':')[1]
                to_save = list(compress(entry, (1,0,1,1,1,0,1,1,0)))
                if transcript_id not in te_map:
                    te_map[transcript_id] = [to_save]
                else:
                    te_map[transcript_id].append(to_save)
            elif parent.startswith('gene:'):
                gene_id = parent.split(':')[1]
                transcript_id = attributes['ID'].split(':')[1]
                if gene_id not in gt_map:
                    gt_map[gene_id] = [transcript_id]
                else:
                    gt_map[gene_id].append(transcript_id)
            else:
                raise Exception('Unexpected #1')


# Write the results
print("Result 1 is being written as a joblib object.")
output_path_1 = os.path.join(output_dir, data_repo_dir, OUTPUT_FILE_1)
joblib.dump(te_map, output_path_1)
print(f"Done: {output_path_1}")

print("Result 2 is being written as a joblib object.")
output_path_2 = os.path.join(output_dir, data_repo_dir, OUTPUT_FILE_2)
joblib.dump(gt_map, output_path_2)
print(f"Done: {output_path_2}")


# End of the script
