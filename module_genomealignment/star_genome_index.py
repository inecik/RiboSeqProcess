"""
Star alignment for paired end sequencing.

"""


import os
import re
import sys
import subprocess
from multiprocessing import cpu_count
from datetime import datetime
from shutil import which


# CONSTANTS
TEMP_DATA_REPO = "star_module"
GENOME_INDEX_DIR_NAME = "genome_index"
ReadLength_minus_1 = 100  # Manual says: In most cases, the default value of 100 will work as well as the ideal value.


# Inputs
temp_dir = sys.argv[1]  # Command line input 4 for the directory for temp files, as described in main.py


# Create and set working directory for genome index creation
temp_dir_module = os.path.join(temp_dir, TEMP_DATA_REPO)
if not os.access(temp_dir_module, os.W_OK) or not os.path.isdir(temp_dir_module):  # Create directory if not exist
    os.mkdir(temp_dir_module)
genome_index_dir = os.path.join(temp_dir_module, GENOME_INDEX_DIR_NAME)
if not os.access(genome_index_dir, os.W_OK) or not os.path.isdir(genome_index_dir):  # Create directory if not exist
    os.mkdir(genome_index_dir)  # For genome index files
os.chdir(temp_dir_module)  # Since everything will be output there for a while


# Download necessary genome files from the server
# Release-96 (GRCh38.p12) is used to be compatible with Mati and Kai's previous works
db_nm_fa = "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fb_nm_gtf = "Homo_sapiens.GRCh38.96.gtf.gz"
db_nm_fa_uncmprsd = re.search(r"(.*)\.fa\.gz$", db_nm_fa).group(1) + ".fa"
fb_nm_gtf_uncmprsd = re.search(r"(.*)\.gtf\.gz$", fb_nm_gtf).group(1) + ".gtf"


# Download genome fasta file
if not os.access(db_nm_fa_uncmprsd, os.R_OK) or not os.path.isfile(db_nm_fa_uncmprsd):
    if not os.access(db_nm_fa, os.R_OK) or not os.path.isfile(db_nm_fa):
        subprocess.run(f"curl -L -R -O ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/{db_nm_fa}", shell=True)
    subprocess.run(f"gzip -d {db_nm_fa}", shell=True)  # Uncompress gz file

# todo: <<<< gff3 kullan >>>>
# Download annotation file
if not os.access(fb_nm_gtf_uncmprsd, os.R_OK) or not os.path.isfile(fb_nm_gtf_uncmprsd):
    if not os.access(fb_nm_gtf, os.R_OK) or not os.path.isfile(fb_nm_gtf):
        subprocess.run(f"curl -L -R -O ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/{fb_nm_gtf}", shell=True)
    subprocess.run(f"gzip -d {fb_nm_gtf}", shell=True)  # Uncompress gz file


# Genome index creation
subprocess.run((
   f"{which('STAR')} "  # Define which star installation to use
   f"--runThreadN {cpu_count()} "  # Define how many core to be used. All cores are now using
   "--runMode genomeGenerate "
   f"--genomeDir {genome_index_dir} "  # Directory to save the files
   f"--genomeFastaFiles {db_nm_fa_uncmprsd} "  # Specifies FASTA file with the genome reference sequences
   f"--sjdbGTFfile {fb_nm_gtf_uncmprsd} "  # Specifies the file with annotated transcripts in the standard GTF format
   f"--sjdbOverhang {ReadLength_minus_1}"  # Specifies the length of the genomic sequence around the annotated junction
), shell=True)


# End of the script
