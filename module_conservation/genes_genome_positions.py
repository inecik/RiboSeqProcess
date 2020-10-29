"""
This script reads GERP scores and single-end SAM file.
Finds out the conservation score of the reads and genomic positions.
It requires 150 GBs of RAM
"""


import os
import subprocess
import joblib


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
GFF3_FILE = "Homo_sapiens.GRCh38.96.chr.gff3"
OUTPUT_DATA_REPO = "gerp_genes"
TEMP_DATA_REPO = "gerp"
OUTPUT_FILE = "genes_genome_positions.joblib"


# Create dir if not exist
data_repo_dir = os.path.join(output_dir, OUTPUT_DATA_REPO)
if not os.access(data_repo_dir, os.W_OK) or not os.path.isdir(data_repo_dir):  # Create directory if not exist
    os.mkdir(data_repo_dir)
# Create dir if not exist
temp_repo_dir = os.path.join(temp_dir, TEMP_DATA_REPO)
if not os.access(temp_repo_dir, os.W_OK) or not os.path.isdir(temp_repo_dir):  # Create directory if not exist
    os.mkdir(temp_repo_dir)


# Create or load SQL database
gff_path = os.path.join(temp_dir, TEMP_DATA_REPO, GFF3_FILE)
if not os.access(gff_path, os.R_OK) or not os.path.isfile(gff_path):
    subprocess.run((f"cd {temp_repo_dir}; curl -L -R -O ftp://ftp.ensembl.org/pub/release-96/gff3/"
                   "homo_sapiens/Homo_sapiens.GRCh38.96.chr.gff3.gz"), shell=True)
    subprocess.run(f"cd {temp_repo_dir}; gzip -d Homo_sapiens.GRCh38.96.chr.gff3.gz", shell=True)


with open(gff_path, "r") as gff_handle:
    gff_raw = gff_handle.readlines()
gff = list()
te_map = dict()  # Transcript genomic element map
gt_map = dict()  # Gene transcript map
for entry in gff_raw:
    if not entry.startswith('#'):
        entry = entry.strip().split('\t')
        attributes = dict([i.split('=') for i in entry[8].split(';')])
        if entry[2] in ["exon", "five_prime_UTR", "three_prime_UTR"]:
            parent = attributes['Parent']
            if parent.startswith('transcript:'):
                transcript_id = parent.split(':')[1]
                if transcript_id not in te_map:
                    te_map[transcript_id] = [[entry[0], entry[2], entry[3], entry[4], entry[6]]]
                else:
                    te_map[transcript_id].append([entry[0], entry[2], entry[3], entry[4], entry[6]])
            else:
                raise Exception('Unexpected #1')
        # Create / add to a dictionary of the mapping
        elif entry[2] == "mRNA":
            parent = attributes['Parent']
            if parent.startswith('gene:'):
                gene_id = parent.split(':')[1]
                transcript_id = attributes['ID'].split(':')[1]
                if gene_id not in gt_map:
                    gt_map[gene_id] = [transcript_id]
                else:
                    gt_map[gene_id].append(transcript_id)
            else:
                raise Exception('Unexpected #2')


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
output_path = os.path.join(output_dir, data_repo_dir, OUTPUT_FILE)
joblib.dump(gene_pos_dict, output_path)


# End of the script
