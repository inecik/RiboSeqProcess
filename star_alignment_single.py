"""
Star alignment for sequencing with one read.

"""


import os
import sys
import subprocess
from multiprocessing import cpu_count
from shutil import which


# CONSTANTS
OUTPUT_DATA_REPO = "star_paired_module"
TEMP_DATA_REPO = "star_paired_module"
GENOME_INDEX_DIR_NAME = "genome_index"
ReadLength_minus_1 = 100  # Manual says: In most cases, the default value of 100 will work as well as the ideal value.


# Inputs
read_file = sys.argv[1]  # Command line input 1 for read 1.
output_dir = sys.argv[2]  # Command line input 2 for the directory for output files, as described in main.py
temp_dir = sys.argv[3]  # Command line input 3 for the directory for temp files, as described in main.py


# Download necessary genome files from the server
# Release-96 (GRCh38.p12) is used to be compatible with Mati and Kai's previous works
db_nm_fa = "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fb_nm_gtf = "Homo_sapiens.GRCh38.96.gtf.gz"
subprocess.run(f"curl -L -R -O ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/{db_nm_fa}", shell=True)
subprocess.run(f"curl -L -R -O ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/{fb_nm_gtf}", shell=True)
subprocess.run(f"gzip -d {db_nm_fa}", shell=True)  # Uncompress gz file
subprocess.run(f"gzip -d {fb_nm_gtf}", shell=True)  # Uncompress gz file


# Create and set working directory for genome index creation
temp_dir_module = os.path.join(temp_dir, TEMP_DATA_REPO)
if not os.access(temp_dir_module, os.W_OK) or not os.path.isdir(temp_dir_module):  # Create directory if not exist
    os.mkdir(temp_dir_module)
os.chdir(temp_dir_module)  # Since everything will be output there for a while
genome_index_dir = os.path.join(temp_dir_module, GENOME_INDEX_DIR_NAME)
if not os.access(genome_index_dir, os.W_OK) or not os.path.isdir(genome_index_dir):  # Create directory if not exist
    os.mkdir(genome_index_dir)  # For genome index files


# Genome index creation
subprocess.run((
    f"{which('star')} "  # Define which star installation to use
    f"--runThreadN {cpu_count()} "  # Define how many core to be used. All cores are now using
    "--runMode genomeGenerate " 
    f"--genomeDir {genome_index_dir} "  # Directory to save the files 
    f"--genomeFastaFiles {db_nm_fa} "  # Specifies FASTA file with the genome reference sequences
    f"--sjdbGTFfile {fb_nm_gtf} "  # Specifies the file with annotated transcripts in the standard GTF format
    f"--sjdbOverhang {ReadLength_minus_1}"  # Specifies the length of the genomic sequence around the annotated junction
), shell=True)


# Create and set working directory for alignment
output_dir_module = os.path.join(output_dir, OUTPUT_DATA_REPO)
if not os.access(output_dir_module, os.W_OK) or not os.path.isdir(output_dir_module):  # Create directory if not exist
    os.mkdir(output_dir_module)
os.chdir(output_dir_module)  # Since everything will be output there


# Actual alignment for sequencing
subprocess.run((
    f"{which('star')} "  # Define which star installation to use
    f"--runThreadN {cpu_count()} "  # Define how many core to be used. All cores are now using
    f"--genomeDir {genome_index_dir} "  # Directory for genome index which has been just created
    f"--readFilesIn {read_file} "
    # All parameters were inherited from Mati-Kai's pipeline.
    "--outFilterMultimapNmax 1 "  
    # "--peOverlapNbasesMin 6 " It is for paired end sequencing
    # "--peOverlapMMp 0.1 " It is for paired end sequencing
    "--outFilterType BySJout "
    "--alignIntronMin 5 "
    f"--outFileNamePrefix {output_dir_module} "
    "--outReadsUnmapped Fastx "
    "--outSAMtype BAM SortedByCoordinate "
    "--outSAMattributes All XS "
    "--quantMode GeneCounts "
    "--twopassMode Basic"
), shell=True)


# End of the script