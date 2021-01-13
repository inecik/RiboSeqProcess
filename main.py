"""
Analysis of Paired End Sequencing of Disome Profiling

The pipeline works only for a pair of reads.
If there is more than o sample or replicates, you have to run the script multiple times.
"""


import os
import sys
import joblib
from archieve.common_functions import *
from commands import *


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# Obtain the root dir for scripts to run properly
scripts_directory = os.path.dirname(__file__)  # Where this package is


# Get inputs for output and temp directories
selection = int(sys.argv[1])
output_dir = sys.argv[2]
temp_dir = sys.argv[3]


# Single-end sequencing
if selection == 1:

    assert len(sys.argv) == 5
    read1 = sys.argv[4]

    # Module 01: Preprocessing
    preprocessing_single(scripts_directory, read1, output_dir, temp_dir)

    # Module 02: Cleanup
    cleanup_input = joblib.load(os.path.join(output_dir, ".01_preprocessing.joblib"))
    cleanup_create_index(scripts_directory, temp_dir)
    cleanup_single(scripts_directory, cleanup_input, output_dir, temp_dir)

    # Module 04: Genome Alignment
    genomealignment_input = joblib.load(os.path.join(output_dir, ".02_cleanup.joblib"))
    genomealignment_create_index(scripts_directory, temp_dir)
    genomealignment_single(scripts_directory, genomealignment_input, output_dir, temp_dir)
    genomealignment_umidedup_single(scripts_directory, output_dir)

    # Module 05: Assignment
    assignment_input = joblib.load(os.path.join(output_dir, ".04_genomealignment.joblib"))
    assignment_ilia(scripts_directory, assignment_input["sam"], output_dir, temp_dir)

# Paired-end sequencing
elif selection in [2, 3]:

    assert len(sys.argv) == 6

    read1 = sys.argv[4]
    read2 = sys.argv[5]

    # Module 01: Preprocessing
    preprocessing_paired(scripts_directory, read1, read2, output_dir, temp_dir)

    # Module 02: Cleanup
    cleanup_inputs = joblib.load(os.path.join(output_dir, ".01_preprocessing.joblib"))
    cleanup_create_index(scripts_directory, temp_dir)
    cleanup_paired(scripts_directory, cleanup_inputs[0], cleanup_inputs[1], output_dir, temp_dir)

    if selection == 2:  # No link pairs

        # Module 04: Genome Alignment
        genomealignment_inputs = joblib.load(os.path.join(output_dir, ".02_cleanup.joblib"))
        genomealignment_create_index(scripts_directory, temp_dir)
        genomealignment_paired(scripts_directory, genomealignment_inputs[0], genomealignment_inputs[1], output_dir, temp_dir)
        genomealignment_umidedup_paired(scripts_directory, output_dir)

    elif selection == 3:  # Link pairs

        # Module 03: Link Pairs
        linkpair_inputs = joblib.load(os.path.join(output_dir, ".02_cleanup.joblib"))
        linkpairs_create_index(scripts_directory, temp_dir)
        linkpairs_alignment(scripts_directory, linkpair_inputs[0], linkpair_inputs[1], output_dir, temp_dir)
        linkpairs_fasta(scripts_directory, output_dir, temp_dir)

        # Module 04: Genome Alignment
        genomealignment_input = joblib.load(os.path.join(output_dir, ".03_linkpairs.joblib"))
        genomealignment_create_index(scripts_directory, temp_dir)
        genomealignment_single(scripts_directory, genomealignment_input, output_dir, temp_dir)
        genomealignment_umidedup_single(scripts_directory, output_dir)

        # Module 05: Assignment
        assignment_input = joblib.load(os.path.join(output_dir, ".04_genomealignment.joblib"))
        assignment_ilia(scripts_directory, assignment_input["sam"], output_dir, temp_dir)


# End of the script
