"""
This script is to compare the results of two bam files on STAR.
It is to test the link pairing module's success.
"""


import os
import sys
import re
import pysam
from pyensembl import EnsemblRelease  #pyensembl install --release 96 --species human

sys.path.append(os.path.abspath(__file__).split("Ribo-seq-Analysis")[0] + "Ribo-seq-Analysis")
from archieve.common_functions import progressBarForTerminal


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# Constants
sam1 = sys.argv[1]  # Single
sam2 = sys.argv[2]  # Paired
try:
    data = EnsemblRelease(96)
except:
    quit(f"Run following command before the code:\npyensembl install --release 96 --species human\n")


# Take the single end sequencing into random access memory
aligned_single = dict()
counter = 0
with pysam.AlignmentFile(sam1, "r") as handle1:  # Open sam files
    sam_iterator = handle1.fetch()  # Get the iterator
    for e in sam_iterator:  # Iterate over the entries
        counter += 1
        qn = re.search(r"^[^_]*", e.query_name).group()  # Remove UMI info if exists
        if qn in aligned_single:
            raise Exception("Multiple entry for the same identifier")
        else:
            aligned_single[qn] = [e.reference_name, e.reference_start]  # Save the entry


# For progress bar
number_of_line_in_sam2 = 0
with pysam.AlignmentFile(sam2, "r") as handle2:  # Open sam files
    sam_iterator = handle2.fetch()  # Get the iterator
    for e in sam_iterator:  # Iterate over the entries
        number_of_line_in_sam2 += 1
# Counters
counter_not_proper_pair = 0
counter_protein_coding_notinsingle = 0
counter_noncoding_notinsingle = 0
reporter = 0
# Open the paired-end file and get only the ones found in single end
aligned_paired_insingle = dict()
with pysam.AlignmentFile(sam2, "r") as handle2:  # Open sam files
    sam_iterator = handle2.fetch()  # Get the iterator
    for e in sam_iterator:  # Iterate over the entries
        reporter += 1
        qn = re.search(r"^[^_]*", e.query_name).group()  # Remove UMI info if exists
        # Check if it is a part of non-proper pair.
        if e.is_unmapped or not e.is_paired or not e.is_proper_pair or e.mate_is_unmapped:
            counter_not_proper_pair += 1
        # Below two 'elif' is to save entries, which are in single-end, to the RAM
        elif qn in aligned_single and qn not in aligned_paired_insingle:
            aligned_paired_insingle[qn] = [e.reference_name, e.reference_start, None]
        elif qn in aligned_single:  # and qn not in aligned_paired_insingle:
            assert aligned_paired_insingle[qn][2] == None, "More than two pair with the same identifier!"
            aligned_paired_insingle[qn][2] = e.reference_start
        else: # Look the ones found in exclusively paired-end
            fepe = data.genes_at_locus(contig=e.reference_name, position=e.reference_start)
            if len(fepe) > 0 and any([i.is_protein_coding for i in fepe]):
                counter_protein_coding_notinsingle += 1
            else:
                counter_noncoding_notinsingle += 1
        if reporter == number_of_line_in_sam2 or reporter % 1000 == 0:
            progressBarForTerminal(reporter, number_of_line_in_sam2)


# Counters
counter_verified_start_position = 0
counter_close_offset = 0
counter_far_apart = 0
counter_both_non_coding = 0
counter_pairedprotein_not_single = 0
counter_singleprotein_not_paired = 0
counter_both_coding = 0
counter_protein_coding_exclusively_single = 0
counter_noncoding_exclusively_single = 0
OFFSET = 1250

for iteration, qn in enumerate(aligned_single): # Iterate through the single end sequencing file.
    se = aligned_single[qn]  # The values corresponding to the query name
    if qn in aligned_paired_insingle:  # Look the ones found in both
        pe = aligned_paired_insingle[qn]  # The values corresponding to the query name

        if pe[0] == se[0] and (se[1] == pe[1] or se[1] == pe[2]): # Check if there is a match
            counter_verified_start_position += 1

        elif None in pe:
            raise Exception("Unexpected! #1")

        elif (pe[0] == se[0] # If start sites are matching for any of pairs and the positions are very close
              and ((-OFFSET < pe[1] - se[1] < OFFSET) or (-OFFSET < pe[2] - se[1] < OFFSET))):
            counter_close_offset += 1

        else:
            gs = data.genes_at_locus(contig=se[0], position=se[1])  # Get the annotation of the position
            gp1 = data.genes_at_locus(contig=pe[0], position=pe[1])  # Do not check the other pair

            counter_far_apart += 1

            if len(gs) <= 0 and len(gp1) <= 0:
                counter_both_non_coding += 1

            elif len(gs) <= 0 and len(gp1) > 0:
                if any([i.is_protein_coding for i in gp1]):
                    counter_pairedprotein_not_single += 1
                else:
                    counter_both_non_coding += 1

            elif len(gs) > 0 and len(gp1) <= 0:
                if any([i.is_protein_coding for i in gs]):
                    counter_singleprotein_not_paired += 1
                else:
                    counter_both_non_coding += 1

            elif len(gs) > 0 and len(gp1) > 0:
                if any([i.is_protein_coding for i in gs]) and any([i.is_protein_coding for i in gp1]):
                    counter_both_coding += 1
                elif any([i.is_protein_coding for i in gs]): # and not any([i.is_protein_coding for i in gp1]):
                    counter_singleprotein_not_paired += 1
                elif any([i.is_protein_coding for i in gp1]):
                    counter_pairedprotein_not_single += 1
                else:
                    counter_both_non_coding += 1

            else:
                raise Exception("Unexpected! #2")

    else: # Look the ones exclusively found in single-end
        gs = data.genes_at_locus(contig=se[0], position=se[1])
        if len(gs) > 0 and any([i.is_protein_coding for i in gs]):
            counter_protein_coding_exclusively_single += 1
        else:
            counter_noncoding_exclusively_single += 1

    if len(aligned_single) - 1 == iteration or iteration % 1000 == 0:
        progressBarForTerminal(iteration, len(aligned_single) - 1)


print("\n")

print(f"Number of entries in single-end sequencing sam: {counter}")
print(f"Number of entries in paired-end sequencing sam: {reporter} ({int(reporter / 2)} pairs)")
print()

print(f"{counter_not_proper_pair} / {reporter} entries in paired-end sam were not in 'proper-pair' or not aligned properly.")
print(f"{len(aligned_paired_insingle)} / {int(reporter / 2)} pairs were in single-end seq sam.")
print(f"{int((counter_protein_coding_notinsingle + counter_noncoding_notinsingle)/2)} / {int(reporter / 2)} pairs were exclusively found in paired-end sam.")
print(f"\t- Among those, {int(counter_protein_coding_notinsingle/2)} pairs were protein coding.")
print(f"\t- Among those, {int(counter_noncoding_notinsingle/2)} pairs were not protein coding.")
print()


print(f"{counter_verified_start_position} / {counter} (%{round(counter_verified_start_position / counter * 100, 2)}) single-end entry has the same start position with one of the pairs (Verified)")
print(f"{len(aligned_single) - len(aligned_paired_insingle)} (%{round((len(aligned_single) - len(aligned_paired_insingle)) / counter * 100, 2)}) reads were found exclusively in single-sam but not paired.")
print(f"\t- Among those, {counter_protein_coding_exclusively_single} reads were protein coding.")
print(f"\t- Among those, {counter_noncoding_exclusively_single} reads were not protein coding.")
print(f"{(counter_far_apart + counter_close_offset)} (%{round((counter_far_apart + counter_close_offset) / counter * 100, 2)}) reads were found in both but could not be verified.")
print(f"\t- Among those, {counter_close_offset} reads were aligned very close (within {OFFSET*2} nt) in single-end compared to paired-end.")
print(f"\t- Among those, {counter_far_apart} reads were aligned significantly differently.")
print(f"\t\t. Both non-coding: {counter_both_non_coding}")
print(f"\t\t. Both protein coding: {counter_both_coding}")
print(f"\t\t. Only single-end protein coding: {counter_singleprotein_not_paired}")
print(f"\t\t. Only paired-end protein coding: {counter_pairedprotein_not_single}")
print()


# End of the script
