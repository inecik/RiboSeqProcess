"""
This script reads GERP scores and single-end SAM file.
Finds out the conservation score of the reads and genomic positions.
It requires 150 GBs of RAM
"""


import os
import sys
import pyBigWig
import subprocess
import pysam
import pandas as pd
import joblib


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# Inputs
sam_path = "/home/kai/KEMALINECIK/out/bams/single.sam"
output_dir = "/home/kai/KEMALINECIK/out/OUTPUT"
temp_dir = "/home/kai/KEMALINECIK/out/TEMP"


# CONSTANTS
GERP_BIGWIG = "gerp_conservation_scores.homo_sapiens.GRCh38.bw"
OUTPUT_DATA_REPO = "gerp_footprints"
TEMP_DATA_REPO = "gerp"
OUTPUT_FILE = "footprint_positions_cons-scores.joblib"


def progressBarForTerminal (iteration, total, prefix ='Progress:', suffix ='', decimals = 1, barLength = 50):
    """
    This function should be called inside of loop, gives the loop's progress.
    :param iteration: It is integer. It is current iteration.
    :param total: It is integer. It is total iteration.
    :param prefix: It is string. It will be placed before progress bar.
    :param suffix: It is string. It will be placed after progress bar.
    :param decimals: It is integer. It is number of decimals in percent complete.
    :param barLength: It is integer. It is character length of bar.
    :return: It is void function. Nothing is returned.
    """
    filledLength = int(round(barLength * iteration / float(total)))
    percents = round(100.00 * (iteration / float(total)), decimals)
    bar = '█' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        sys.stdout.write('\n')
        sys.stdout.flush()


# Create dir if not exist
data_repo_dir = os.path.join(output_dir, OUTPUT_DATA_REPO)
if not os.access(data_repo_dir, os.W_OK) or not os.path.isdir(data_repo_dir):  # Create directory if not exist
    os.mkdir(data_repo_dir)
# Create dir if not exist
temp_repo_dir = os.path.join(temp_dir, TEMP_DATA_REPO)
if not os.access(temp_repo_dir, os.W_OK) or not os.path.isdir(temp_repo_dir):  # Create directory if not exist
    os.mkdir(temp_repo_dir)


# Download genome conservation file
gerp_bigwig_path = os.path.join(temp_repo_dir, GERP_BIGWIG)
if not os.access(gerp_bigwig_path, os.R_OK) or not os.path.isfile(gerp_bigwig_path):  # Download if not exists
    print(gerp_bigwig_path)
    print("BigWig file is downloading.")
    subprocess.run((f"cd {temp_repo_dir}; curl -L -R -O ftp://ftp.ensembl.org/pub/release-96/compara/"
                   "conservation_scores/88_mammals.gerp_conservation_score/"
                   "gerp_conservation_scores.homo_sapiens.GRCh38.bw"), shell=True)


# Take all conservation file into memory
cons = dict()
with pyBigWig.open(gerp_bigwig_path, "r") as bw, pysam.AlignmentFile(sam_path, "r") as sam:
    chr_mem = [i for i in sam.references if i in list(bw.chroms())]
    print("Reading conservation scores of chromosome:")
    for chr in chr_mem:
        print(chr)
        cons[chr] = bw.values(chr, 0, bw.chroms()[chr])
print("Score reading completed")


# For progress bar
number_of_line_in_sam = 0
with pysam.AlignmentFile(sam_path, "r") as sam:  # Open sam files
    sam_iterator = sam.fetch()  # Get the iterator
    for e in sam_iterator:  # Iterate over the entries
        number_of_line_in_sam += 1


# Fetch the score and crete a dataframe to store the data
df_footprints_out = pd.DataFrame(columns=["query_name", "chromosome", "positions", "scores"])
reporter = 0
with pyBigWig.open(gerp_bigwig_path, "r") as bw, pysam.AlignmentFile(sam_path, "r") as sam:
    sam_iterator = sam.fetch()  # Get the iterator
    for e in sam_iterator:  # Iterate over the entries
        reporter += 1
        positions = e.get_reference_positions()
        scores = [cons[e.reference_name][i] for i in positions]
        df_footprints_out.append({
            "query_name": e.query_name,
            "chromosome": e.reference_name,
            "positions": positions,
            "scores": scores,
        }, ignore_index=True)
        if reporter == number_of_line_in_sam or reporter % 10000 == 0:
            progressBarForTerminal(reporter, number_of_line_in_sam)


# Write the results
output_path = os.path.join(data_repo_dir, OUTPUT_FILE)
joblib.dump(df_footprints_out, output_path)


# End of the script
