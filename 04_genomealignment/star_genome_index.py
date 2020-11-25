"""
Star alignment for paired end sequencing.

"""


import os
import sys
import subprocess
from multiprocessing import cpu_count
from shutil import which

sys.path.append(os.path.abspath(__file__).split("Ribo-seq-Analysis")[0] + "Ribo-seq-Analysis")
from supplementary.common_functions import bcolors as c


# CONSTANTS
TEMP_DATA_REPO = "star_module"
GENOME_INDEX_DIR_NAME = "genome_index"
ReadLength_minus_1 = 100  # Manual says: In most cases, the default value of 100 will work as well as the ideal value.


# Inputs
temp_dir = sys.argv[1]  # Command line input 4 for the directory for temp files, as described in main.py

# Check the genome files are there or not.
print("Checking genome indexing integrity. It is just based on file names.")
temp_dir_module = os.path.join(temp_dir, TEMP_DATA_REPO)
genome_index_dir = os.path.join(temp_dir_module, GENOME_INDEX_DIR_NAME)
if os.access(genome_index_dir, os.R_OK) or os.path.isdir(genome_index_dir):
    files_reference = ['Log.out', 'exonGeTrInfo.tab', 'geneInfo.tab', 'transcriptInfo.tab', 'exonInfo.tab', 'Genome',
                       'sjdbList.fromGTF.out.tab', 'chrName.txt', 'chrStart.txt', 'chrLength.txt', 'chrNameLength.txt',
                       'sjdbInfo.txt', 'sjdbList.out.tab', 'genomeParameters.txt', 'SA', 'SAindex']
    print(genome_index_dir)
    files_in_dir = os.listdir(genome_index_dir)
    if sorted(files_reference) == sorted(files_in_dir):
        print(f"{c.WARNING}Genome indexing files are found. Genome indexing is skipped.{c.ENDC}")
        sys.exit()
    else:
        print(f"{c.WARNING}Genome indexing files are not found.{c.ENDC}")

# Create and set working directory for genome index creation
if not os.access(temp_dir_module, os.W_OK) or not os.path.isdir(temp_dir_module):  # Create directory if not exist
    os.mkdir(temp_dir_module)
    print("Temp directory created.")
if not os.access(genome_index_dir, os.W_OK) or not os.path.isdir(genome_index_dir):  # Create directory if not exist
    os.mkdir(genome_index_dir)  # For genome index files
    print("Genome indexing directory created under temp.")
os.chdir(temp_dir_module)  # Since everything will be output there for a while


# Download necessary genome files from the server
db_nm_fa = "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fb_nm_gtf = "Homo_sapiens.GRCh38.100.gtf.gz"
db_nm_fa_uncmprsd = os.path.splitext(db_nm_fa)[0]
fb_nm_gtf_uncmprsd = os.path.splitext(fb_nm_gtf)[0]


# Download genome fasta file
print("Human genome fasta file is being downloading from the server.")
if not os.access(db_nm_fa_uncmprsd, os.R_OK) or not os.path.isfile(db_nm_fa_uncmprsd):
    if not os.access(db_nm_fa, os.R_OK) or not os.path.isfile(db_nm_fa):
        subprocess.run(f"curl -L -R -O ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/{db_nm_fa}", shell=True)
    subprocess.run(f"gzip -d {db_nm_fa}", shell=True)  # Uncompress gz file


# Download annotation file
print("Human genome gtf file is being downloading from the server.")
if not os.access(fb_nm_gtf_uncmprsd, os.R_OK) or not os.path.isfile(fb_nm_gtf_uncmprsd):
    if not os.access(fb_nm_gtf, os.R_OK) or not os.path.isfile(fb_nm_gtf):
        subprocess.run(f"curl -L -R -O ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/{fb_nm_gtf}", shell=True)
    subprocess.run(f"gzip -d {fb_nm_gtf}", shell=True)  # Uncompress gz file


# Genome index creation
print("Genome indexing is started.")
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
