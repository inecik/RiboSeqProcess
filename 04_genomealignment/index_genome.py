"""
Star alignment for paired end sequencing.

"""


import os
import sys
import subprocess
from multiprocessing import cpu_count
from shutil import which
import joblib

repository_name = "RiboSeqProcess"
sys.path.append(os.path.abspath(__file__).split(repository_name)[0] + repository_name)
from archieve.common_functions import *


# CONSTANTS
TEMP_DATA_REPO = "04_genomealignment"
GENOME_INDEX_DIR_NAME = "index_genome"
ReadLength_minus_1 = 100  # Manual says: "In most cases, the default value of 100 will work as well as the ideal value."


# Inputs
temp_dir = sys.argv[1]  # as described in main.py
temp_dir_module = create_dir(temp_dir, TEMP_DATA_REPO)
genome_index_dir = create_dir(temp_dir_module, GENOME_INDEX_DIR_NAME)
metadata_index_path = os.path.join(temp_dir_module, ".04_genomealignment_metadata.joblib")
os.chdir(temp_dir_module)  # Since everything will be output there for a while


try:
    # Check if the index file
    assert os.path.isfile(metadata_index_path) and os.access(metadata_index_path, os.R_OK)
    metadata_index_previously = joblib.load(metadata_index_path)
    metadata_index = get_files_metadata(genome_index_dir)
    assert metadata_index_previously == metadata_index
    print("Index files for genome alignment are already present.")

except:
    print("Index files for genome alignment are now being created.")

    # Download necessary genome files from the server
    fa = "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    gtf = "Homo_sapiens.GRCh38.100.gtf.gz"
    fa_uncompressed = os.path.splitext(fa)[0]
    gtf_uncompressed = os.path.splitext(gtf)[0]

    # Download genome fasta file
    print("Human genome fasta file is being downloading from the server.")
    if not os.access(fa_uncompressed, os.R_OK) or not os.path.isfile(fa_uncompressed):
        if not os.access(fa, os.R_OK) or not os.path.isfile(fa):
            subprocess.run(f"curl -L -R -O ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/{fa}", shell=True)
        subprocess.run(f"gzip -d {fa}", shell=True)  # Uncompress gz file

    # Download annotation file
    print("Human genome gtf file is being downloading from the server.")
    if not os.access(gtf_uncompressed, os.R_OK) or not os.path.isfile(gtf_uncompressed):
        if not os.access(gtf, os.R_OK) or not os.path.isfile(gtf):
            subprocess.run(f"curl -L -R -O ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/{gtf}", shell=True)
        subprocess.run(f"gzip -d {gtf}", shell=True)  # Uncompress gz file

    # Genome index creation
    print("Genome indexing is started.")
    subprocess.run((
       f"{which('STAR')} "  # Define which star installation to use
       f"--runThreadN {cpu_count()} "  # Define how many core to be used. All cores are now using
       "--runMode genomeGenerate "
       f"--genomeDir {genome_index_dir} "  # Directory to save the files
       f"--genomeFastaFiles {fa_uncompressed} "  # Specifies FASTA file with the genome reference sequences
       f"--sjdbGTFfile {gtf_uncompressed} "  # Specifies the file with annotated transcripts in the standard GTF format
       f"--sjdbOverhang {ReadLength_minus_1}"  # Specifies the length of the genomic sequence around the annotated junction
    ), shell=True)

    # Write the file info
    metadata_index = get_files_metadata(genome_index_dir)
    joblib.dump(metadata_index, metadata_index_path)

# End of the script
