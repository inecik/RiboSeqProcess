"""
Ilia's julia script 02_assign_count_streaming_11.10.19.jl fails when there are multiple entries
for a given gene name. This script renames gene names if there are more than one occurences.
"""


import os
import sys


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


# Inputs
gff_path = sys.argv[1]
output_path = os.path.splitext(gff_path)[0] + "_renamed_duplicate_gene_names.gff3"


print("GFF3 is imported to the random access memory.")
with open(gff_path, "r") as gff_handle:
    gff_raw = gff_handle.readlines()

print("Gene names are read.")
n_map = dict()
for entry in gff_raw:
    if not entry.startswith('#'):
        entry = entry.strip().split('\t')
        attributes = dict([i.split('=') for i in entry[8].split(';')])
        if entry[2] == "gene":
            gene_name = attributes['Name']
            gene_id = attributes['ID'].split(':')[1]
            if gene_name not in n_map:
                n_map[gene_name] = [gene_id]
            else:
                n_map[gene_name].append(gene_id)


print("Duplicated gene names are detected.")
duplicated_gene_name = dict()
for i in n_map:
    if len(n_map[i]) > 1:
        duplicated_gene_name[i] = 0


print("Output is written by changing the names of duplicated gene names.")
with open(output_path, "w") as output_handle:
    for entry in gff_raw:
        if entry.startswith('#'):
            output_handle.write(entry)
        else:
            entry_split = entry.strip().split('\t')
            attributes = dict([i.split('=') for i in entry_split[8].split(';')])
            if entry_split[2] == "gene" and attributes['Name'] in duplicated_gene_name:
                the_name = attributes['Name']
                duplicated_gene_name[the_name] += 1
                new_name = the_name + f"_duplicated_{duplicated_gene_name[the_name]}"
                new_entry = entry.replace(f";Name={the_name};", f";Name={new_name};")
                output_handle.write(new_entry)
            else:
                output_handle.write(entry)


# End of the script
