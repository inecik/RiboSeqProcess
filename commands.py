"""

"""


import os
import subprocess
import sys
from shutil import which
from archieve.common_functions import bcolors


def preprocessing_paired(scripts_directory, read1, read2, output_dir, temp_dir):
    
    print(f"{bcolors.HEADER}Preprocessing Module.{bcolors.ENDC}")

    spr = subprocess.run((
        f"{which('python3')} "  # Define which python installation to use
        f"{os.path.join(scripts_directory, '01_preprocessing/cutadapt_umitools_paired.py')} "  # Relative script dir
        f"{read1} "  # sys.argv[1]
        f"{read2} "  # sys.argv[2]
        f"{output_dir} "
        f"{temp_dir}"
    ), shell=True)
    
    if spr.returncode != 0:
        sys.exit(f"{bcolors.FAIL}Error in cutadapt_umitools_paired.py Exiting.{bcolors.ENDC}")
    
    print(f"{bcolors.OKCYAN}Preprocessing Module was successfully completed.\n{bcolors.ENDC}")
    

def preprocessing_single(scripts_directory, read1, output_dir, temp_dir):

    print(f"{bcolors.HEADER}Preprocessing Module.{bcolors.ENDC}")

    spr = subprocess.run((
        f"{which('python3')} "  # Define which python installation to use
        f"{os.path.join(scripts_directory, '01_preprocessing/cutadapt_umitools_single.py')} "  # Relative script dir
        f"{read1} "  # sys.argv[1]
        f"{output_dir} "
        f"{temp_dir}"
    ), shell=True)

    if spr.returncode != 0:
        sys.exit(f"{bcolors.FAIL}Error in cutadapt_umitools_single.py: Exiting.{bcolors.ENDC}")

    print(f"{bcolors.OKCYAN}Preprocessing Module was successfully completed.\n{bcolors.ENDC}")


def cleanup_create_index(scripts_directory, temp_dir):
    
    print(f"{bcolors.HEADER}Clean-up Module indexing.{bcolors.ENDC}")
    
    spr = subprocess.run((
        f"{which('python3')} "  # Define which python installation to use
        f"{os.path.join(scripts_directory, '02_cleanup/index_ncrna.py')} "  # Relative script dir
        f"{temp_dir} "  # sys.argv[1]
    ), shell=True)
    
    if spr.returncode != 0: 
        sys.exit(f"{bcolors.FAIL}Error in database_rnaremove_bowtie2.py: Exiting.{bcolors.ENDC}")
    
    print(f"{bcolors.OKCYAN}Genome indexing for rRNA removal was successfully completed.\n{bcolors.ENDC}")


def cleanup_paired(scripts_directory, read1, read2, output_dir, temp_dir):

    print(f"{bcolors.HEADER}Clean-up Module.{bcolors.ENDC}")
    
    spr = subprocess.run((
        f"{which('python3')} "  # Define which python installation to use
        f"{os.path.join(scripts_directory, '02_cleanup/rnaremove_paired.py')} " 
        f"{read1} "  # sys.argv[1] Read 1
        f"{read2} "  # sys.argv[2] Read 2
        f"{output_dir} "  # sys.argv[3]
        f"{temp_dir}"  # sys.argv[3]
    ), shell=True)
    
    if spr.returncode != 0: 
        sys.exit(f"{bcolors.FAIL}Error in bowtie2_rnaremove.py: Exiting.{bcolors.ENDC}")
    
    print(f"{bcolors.OKCYAN}rRNA removal was successfully completed.\n{bcolors.ENDC}")


def cleanup_single(scripts_directory, read1, output_dir, temp_dir):

    print(f"{bcolors.HEADER}Clean-up Module.{bcolors.ENDC}")

    spr = subprocess.run((
        f"{which('python3')} "  # Define which python installation to use
        f"{os.path.join(scripts_directory, '02_cleanup/rnaremove_single.py')} "
        f"{read1} "  # sys.argv[1] Read 1
        f"{output_dir} "  # sys.argv[3]
        f"{temp_dir}"  # sys.argv[3]
    ), shell=True)

    if spr.returncode != 0:
        sys.exit(f"{bcolors.FAIL}Error in bowtie2_rnaremove.py: Exiting.{bcolors.ENDC}")

    print(f"{bcolors.OKCYAN}rRNA removal was successfully completed.\n{bcolors.ENDC}")


def linkpairs_create_index(scripts_directory, temp_dir):
    
    print(f"{bcolors.HEADER}Link Pairing Indexing.{bcolors.ENDC}")

    spr = subprocess.run((
        f"{which('python3')} "  # Define which python installation to use
        f"{os.path.join(scripts_directory, '03_linkpairs/index_transcriptome.py')} "  # Relative script dir
        f"{temp_dir} "  # sys.argv[1]
    ), shell=True)
    
    if spr.returncode != 0: 
        sys.exit(f"{bcolors.FAIL}Error in database_transcriptome_bowtie2.py: Exiting.{bcolors.ENDC}")
    
    print(f"{bcolors.OKCYAN}Genome indexing for linking pairs was successfully completed.\n{bcolors.ENDC}")


def linkpairs_alignment(scripts_directory, read1, read2, output_dir, temp_dir):

    print(f"{bcolors.HEADER}Link Pairing Alignment.{bcolors.ENDC}")
    
    spr = subprocess.run((
        f"{which('python3')} "  # Define which python installation to use
        f"{os.path.join(scripts_directory, '03_linkpairs/prealignment.py')} "
        f"{read1} "  # sys.argv[1] 
        f"{read2} "  # sys.argv[2]
        f"{output_dir} "  # sys.argv[3]
        f"{temp_dir}"  # sys.argv[3]
    ), shell=True)
    
    if spr.returncode != 0: 
        sys.exit(f"{bcolors.FAIL}Error in bowtie2_prealignment.py: Exiting.{bcolors.ENDC}")
    
    print(f"{bcolors.OKCYAN}Alignment for linking pairs was successfully completed.\n{bcolors.ENDC}")


def linkpairs_fasta(scripts_directory, output_dir, temp_dir):

    print(f"{bcolors.HEADER}Link Pairing Create Fasta.{bcolors.ENDC}")

    spr = subprocess.run((
        f"{which('python3')} "  # Define which python installation to use
        f"{os.path.join(scripts_directory, '03_linkpairs/sam_processor.py')} " 
        f"{output_dir} "  # sys.argv[3]
        f"{temp_dir}"  # sys.argv[3]
    ), shell=True)
    
    if spr.returncode != 0: 
        sys.exit(f"{bcolors.FAIL}Error in sam_processor.py: Exiting.{bcolors.ENDC}")
    
    print(f"{bcolors.OKCYAN}Creation of fasta with linked pairs was successfully completed.\n{bcolors.ENDC}")


def genomealignment_create_index(scripts_directory, temp_dir):
    
    print(f"{bcolors.HEADER}Genome Alignment Index.{bcolors.ENDC}")
    
    spr = subprocess.run((
        f"{which('python3')} "  # Define which python installation to use
        f"{os.path.join(scripts_directory, '04_genomealignment/index_genome.py')} "  # Relative script dir
        f"{temp_dir}"  # sys.argv[1]
    ), shell=True)
    
    if spr.returncode != 0:
        sys.exit(f"{bcolors.FAIL}Error in star_genome_index.py: Exiting.{bcolors.ENDC}")
    
    print(f"{bcolors.OKCYAN}Human genome indexing was successfully completed.\n{bcolors.ENDC}")


def genomealignment_single(scripts_directory, read1, output_dir, temp_dir):
    
    print(f"{bcolors.HEADER}Genome Alignment Module.{bcolors.ENDC}")
    
    spr = subprocess.run((
        f"{which('python3')} "  # Define which python installation to use
        f"{os.path.join(scripts_directory, '04_genomealignment/genome_alignment_single.py')} "
        f"{read1} "  # sys.argv[1]
        f"{output_dir} "  # sys.argv[2]
        f"{temp_dir}"  # sys.argv[3]
    ), shell=True)
    
    if spr.returncode != 0:
        sys.exit(f"{bcolors.FAIL}Error in star_alignment_single.py: Exiting.{bcolors.ENDC}")
    
    print(f"{bcolors.OKCYAN}Human genome alignment was successfully completed.\n{bcolors.ENDC}")


def genomealignment_paired(scripts_directory, read1, read2, output_dir, temp_dir):
    pass  # todo


def genomealignment_umidedup_single(scripts_directory, output_dir):

    print(f"{bcolors.HEADER}Genome Alignment UMI deduplication.{bcolors.ENDC}")

    spr = subprocess.run((
        f"{which('python3')} "  # Define which python installation to use
        f"{os.path.join(scripts_directory, '04_genomealignment/umitools_dedup_single.py')} "  # Relative script dir
        f"{output_dir}"
    ), shell=True)
    
    if spr.returncode != 0:
        sys.exit(f"{bcolors.FAIL}Error in umitools_dedup.py: Exiting.{bcolors.ENDC}")

    print(f"{bcolors.OKCYAN}UMI deduplication was successfully completed.\n{bcolors.ENDC}")


def genomealignment_umidedup_paired(scripts_directory, output_dir):
    pass  # todo

def assignment_ilia(scripts_directory, sam_file, output_dir, temp_dir):

    print(f"{bcolors.HEADER}Ilia's Julia script for assignment.{bcolors.ENDC}")

    spr = subprocess.run((
        f"{which('python3')} "  # Define which python installation to use
        f"{os.path.join(scripts_directory, '05_julia_assignment/ilias_assignment_count_umidedup_streaming.py')} "  # Relative script dir
        f"{sam_file} "
        f"{output_dir} "
        f"{temp_dir}"
    ), shell=True)

    if spr.returncode != 0:
        sys.exit(f"{bcolors.FAIL}Error in ilias_assignment_count_umidedup_streaming.py: Exiting.{bcolors.ENDC}")

    print(f"{bcolors.OKCYAN}Ilia's assignment script was successfully completed.\n{bcolors.ENDC}")


# End of the script
