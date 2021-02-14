"""

"""

import sys
import os
import shutil
import subprocess
from Bio import SeqIO
import joblib
import multiprocessing
import pysam
import re
import argparse
import datetime


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.2"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


def main():

    a = argument_parser()

    job_list = JobList(a.filepath)
    job_list.confirm_job_list()
    print(f"{Col.H}Operation started.\n{Col.E}")
    controller = Controller(a.identifier, a.organism, a.ensembl_release, a.temp, a.output, job_list.jobs)
    joblib.dump(controller, os.path.join(a.output, "pipeline_controller.joblib"))
    print(f"{Col.H}Operation successfully ended.\n{Col.E}")

def argument_parser():

    parser = argparse.ArgumentParser(prog="pipeline.py",
                                     description='The app description')  # todo

    parser.add_argument("-r", type=str, dest="identifier",
                        help="The identifier for this run, which will be the directory name for outputs.",
                        required=False, default=datetime.datetime.now().strftime("Run_on_%Y_%m_%d_at_%H_%M_%S"))

    parser.add_argument("-a", type=str, dest="organism",
                        help="Organism of interest for the analysis.",
                        required=False, default="homo_sapiens",
                        choices=["homo_sapiens", "mus_musculus"],  # todo
                        )

    parser.add_argument("-e", type=int, dest="ensembl_release",
                        help="Ensembl version to be used.",
                        required=False, default=102)

    parser.add_argument("-f", type=str, required=True, dest="filepath",
                        help="File path for the task list.")

    parser.add_argument("-t", type=str, required=True, dest="temp",
                        help="Directory to be used for temporary files.")

    parser.add_argument("-o", type=str, required=True, dest="output",
                        help="Directory to be used for output files.")

    return parser.parse_args()


class JobList:

    process_map = {1: "01_Preprocessing",
                   2: "02_Cleanup",
                   3: "03_LinkingPairs",
                   4: "04_GenomeAlignment",
                   5: "05_JuliaAssignment"}

    default_adapter_single = "ATCGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    default_adapter1_paired = "ATCGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    default_adapter2_paired = None
    default_pattern_single = "^(?P<umi_1>.{2}).*(?P<umi_2>.{5})$"
    default_pattern1_paired = "^(?P<umi_1>.{2}).*"
    default_pattern2_paired = ".*(?P<umi_2>.{5})(?P<discard_2>.{5})$"
    default_processes_single = [1, 2, 4, 5]
    default_processes_paired_linking = [1, 2, 3, 4, 5]
    default_processes_paired = [1, 2, 4]
    proper_input_extensions = ["fastq.gz"]

    def __init__(self, user_task_file):
        self.user_task_file = user_task_file
        self.jobs = self.read_job_list(self.user_task_file)

    @staticmethod
    def read_job_list(file_path):

        result = dict()

        with open(file_path, "r") as handle:

            for entry in [r.strip() for r in handle.read().split(">") if r.strip()]:

                if entry.startswith("#"):
                    continue

                line = [r.strip() for r in entry.split(os.linesep) if r.strip()]
                print_error_entry = Col.B + '\n'.join(line) + Col.E

                if len(line[0].split()) != 1:
                    raise ValueError(f"Improper key name:\n{print_error_entry}")
                if line[1] not in ["single", "paired_linking", "paired"]:
                    raise ValueError(f"Improper sequencing method:\n{print_error_entry}")
                if not check_file(line[2], JobList.proper_input_extensions):
                    raise ValueError(f"Improper file or file extension method:\n{print_error_entry}")
                if not os.path.isabs(line[2]):
                    raise ValueError(f"File path has to be absolute:\n{print_error_entry}")
                if line[1].startswith("paired") and not check_file(line[3], JobList.proper_input_extensions):
                    raise ValueError(f"Improper file or file extension method:\n{print_error_entry}")
                if line[1].startswith("paired") and not os.path.isabs(line[3]):
                    raise ValueError(f"File path has to be absolute:\n{print_error_entry}")

                if line[1] == "single":
                    temp_dict = {"sequencing_method": line[1],
                                 "input_fastq": line[2],
                                 "adapter": JobList.default_adapter_single,
                                 "pattern_umi": JobList.default_pattern_single,
                                 "processes": JobList.default_processes_single}
                else:
                    process_temp = JobList.default_processes_paired if line[1] == "paired" \
                        else JobList.default_processes_paired_linking
                    temp_dict = {"sequencing_method": line[1],
                                 "input_fastq": [line[2], line[3]],
                                 "adapter1": JobList.default_adapter1_paired,
                                 "adapter2": JobList.default_adapter2_paired,
                                 "pattern_umi1": JobList.default_pattern1_paired,
                                 "pattern_umi2": JobList.default_pattern2_paired,
                                 "processes": process_temp}

                non_default_entries = line[3:] if line[1] == "single" else line[4:]
                if line[1] == "single":
                    possible_settings = ["adapter", "pattern_umi", "processes"]
                else:
                    possible_settings = ["adapter1", "adapter2", "pattern_umi1", "pattern_umi2" "processes"]
                for nde in non_default_entries:
                    des = [r.strip() for r in nde.split("=") if r.strip()]
                    if len(des) != 2 or des[0] not in possible_settings:
                        print(f"{Col.W}Manual setting ignored for {line[0]}:\n{des}{Col.E}")
                    elif des[0] == "processes":
                        temp_dict[des[0]] = sorted([int(r.strip()) for r in des[1].split(',') if r.strip()])
                    else:
                        temp_dict[des[0]] = des[1] if des[1] != "None" else None

                assert line[0] not in result, f"Duplicated key for {line[0]}"
                result[line[0]] = temp_dict

        return result

    def confirm_job_list(self):
        print(f"{Col.H}Please confirm the jobs.\n"
              f"Press enter to confirm, type anything to stop.\n"
              f"There are {len(self.jobs)} jobs:{Col.E}\n")
        for job_id in self.jobs:
            the_job = self.jobs[job_id]
            processes_pre = ", ".join([str(i) for i in the_job['processes']])

            if the_job["sequencing_method"] == "single":

                adapter_temp = the_job["adapter"] if the_job["adapter"] != JobList.default_adapter_single \
                    else f"Default (\"{the_job['adapter']}\")"
                pattern_temp = the_job["pattern_umi"] if the_job["pattern_umi"] != JobList.default_pattern_single \
                    else f"Default (\"{the_job['pattern_umi']}\")"
                processes_temp = processes_pre if the_job["processes"] != JobList.default_processes_single \
                    else f"Default (\"{processes_pre}\")"

                print_string = (
                    f"{Col.H}> {job_id}{Col.C}{os.linesep}"
                    f"Sequencing Method : {the_job['sequencing_method']}{os.linesep}"
                    f"Read              : {the_job['input_fastq']}{os.linesep}"
                    f"Adapter           : {adapter_temp}{os.linesep}"
                    f"Pattern UMI       : {pattern_temp}{os.linesep}"
                    f"Processes         : {processes_temp}{Col.E}{os.linesep}"
                )

            else:

                adapter1_temp = the_job["adapter1"] if the_job["adapter1"] != JobList.default_adapter1_paired \
                    else f"Default (\"{the_job['adapter1']}\")"
                adapter2_temp = the_job["adapter2"] if the_job["adapter2"] != JobList.default_adapter2_paired \
                    else f"Default (\"{the_job['adapter2']}\")"
                pattern1_temp = the_job["pattern_umi1"] if the_job["pattern_umi1"] != JobList.default_pattern1_paired \
                    else f"Default (\"{the_job['pattern_umi1']}\")"
                pattern2_temp = the_job["pattern_umi2"] if the_job["pattern_umi2"] != JobList.default_pattern2_paired \
                    else f"Default (\"{the_job['pattern_umi2']}\")"

                processes_method = JobList.default_processes_paired if the_job["sequencing_method"] == "paired" \
                    else JobList.default_processes_paired_linking
                processes_temp = processes_pre if the_job["processes"] != processes_method \
                    else f"Default (\"{processes_pre}\")"

                print_string = (
                    f"{Col.H}> {job_id}{Col.C}{os.linesep}"
                    f"Sequencing Method : {the_job['sequencing_method']}{os.linesep}"
                    f"Read 1            : {the_job['input_fastq'][0]}{os.linesep}"
                    f"Read 2            : {the_job['input_fastq'][1]}{os.linesep}"
                    f"Adapter 1         : {adapter1_temp}{os.linesep}"
                    f"Adapter 2         : {adapter2_temp}{os.linesep}"
                    f"Pattern UMI 1     : {pattern1_temp}{os.linesep}"
                    f"Pattern UMI 2     : {pattern2_temp}{os.linesep}"
                    f"Processes         : {processes_temp}{Col.E}{os.linesep}"
                )

            confirm = input(print_string)
            if confirm != "":
                print(f"{Col.F}Process terminated.{Col.E}")
                sys.exit(1)
            else:
                print(f"{Col.W}Confirmed.\n{Col.E}")


class Controller:

    def __init__(self, run_identifier, organism, ensembl_release, temp_repo_dir, data_repo_dir, jobs):

        check_directory([temp_repo_dir, data_repo_dir])
        check_exist_package(["cutadapt", "umi_tools", "bowtie2-build", "STAR", "bowtie2", "/usr/bin/samtools"])

        self.temp_repo_dir = temp_repo_dir
        self.data_repo_dir = data_repo_dir
        self.organism = organism
        self.ensembl_release = ensembl_release
        self.jobs = jobs
        self.run_identifier = run_identifier
        self.org_db = OrganismDatabase(self.organism, self.ensembl_release, self.temp_repo_dir)
        self.cpu = multiprocessing.cpu_count()

        self.julia_path = os.path.join(os.path.dirname(__file__), "julia_assignment.jl")
        self.create_output_tree()
        self.index_directories = self.create_index()

        self.preprocessing()
        self.cleanup()
        self.linking_pairs()
        self.genome_alignment()
        self.julia_assignment()

    def julia_assignment(self):
        for job_id in self.jobs:
            if 5 in self.jobs[job_id]["processes"]:
                self.genome_alignment_one_job(job_id)
                # No output path

    def julia_assignment_one_job(self, job_id):

        jobbing = self.jobs[job_id]
        job_dir = jobbing["processes_dirs"][5]  # 5'i açıkla
        input_sam = self.jobs[job_id]["process_genome_alignment"]
        assert JobList.process_map[4].endswith("JuliaAssignment")

        subprocess.run((
            f"cd {job_dir}; "  # Change the directory to the directory
            f"{shutil.which('julia')} {self.julia_path} "  # Which Julia installation to use and the script
            f"-g {self.julia_assignment_gff3_correct()} "  # Gff3 file. Removed of duplicated gene names
            "-a 3 "  # Assignment from 3'
            # "-u "  # Inherited from Mati. Removed because umi-tool deduplication is already done.
            f"-o {job_dir} {input_sam}"  # Output file & Input file
        ), shell=True)

        # No output path

    def julia_assignment_gff3_correct(self):
        gff_path = self.org_db.get_db("gff3")
        output_path = os.path.splitext(gff_path)[0] + "_renamed_duplicate_gene_names.gff3"

        if not os.access(output_path, os.R_OK) or not os.path.isfile(output_path):
            with open(gff_path, "r") as gff_handle:  # "GFF3 is imported to the random access memory."
                gff_raw = gff_handle.readlines()

            n_map = dict()  # Gene names are read.
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

            duplicated_gene_name = dict()  # Duplicated gene names are detected.
            for i in n_map:
                if len(n_map[i]) > 1:
                    duplicated_gene_name[i] = 0

            # Output is written by changing the names of duplicated gene names.
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
        return output_path

    def genome_alignment(self):
        for job_id in self.jobs:
            if 4 in self.jobs[job_id]["processes"]:  # todo: önceki step yoksa olamıyor olsun
                self.genome_alignment_one_job(job_id)
            else:
                self.jobs[job_id]["process_genome_alignment"] = self.jobs[job_id]["process_linking_pairs"]

    def genome_alignment_one_job(self, job_id):

        jobbing = self.jobs[job_id]
        job_dir = jobbing["processes_dirs"][4]  # 4'i açıkla
        input_fastq_fasta = self.jobs[job_id]["process_linking_pairs"]
        assert JobList.process_map[4].endswith("GenomeAlignment")

        if jobbing["sequencing_method"] in ["paired"]:  # note it down, not paired-linking

            print(f"Genome alignment for {job_id} is now running.")
            prefix = "genome_alignment_"
            subprocess.run((
                f"cd {job_dir}; "  # Change the directory to the index directory
                f"{shutil.which('STAR')} "  # Define which star installation to use
                f"--runThreadN {self.cpu - 2} "  # Define how many core to be used. All cores are now using
                f"--genomeDir {self.index_directories['index_dna']} "
                f"--readFilesIn {input_fastq_fasta[0]} {input_fastq_fasta[1]} "
                # All parameters were inherited from Mati-Kai's pipeline.
                "--outFilterMultimapNmax 1 "
                "--peOverlapNbasesMin 6 "
                "--peOverlapMMp 0.1 "
                "--outFilterType BySJout "
                "--alignIntronMin 5 "
                f"--outFileNamePrefix {os.path.join(job_dir, prefix)} "
                "--outReadsUnmapped Fastx "
                "--outSAMtype BAM SortedByCoordinate "
                "--outSAMattributes All XS "
                "--quantMode GeneCounts "
                "--twopassMode Basic "
                "> report_genome_alignment.log"
            ), shell=True)
            raw_bam = os.path.join(job_dir, prefix + "Aligned.sortedByCoord.out.bam")

            after_sort = os.path.splitext(raw_bam)[0] + '_sorted.bam'
            subprocess.run((  # Sort bam file and Index
                f"cd {job_dir}; "
                f"{shutil.which('/usr/bin/samtools')} sort {raw_bam} -o {after_sort}; " 
                f"{shutil.which('/usr/bin/samtools')} index {after_sort}"
            ), shell=True)

            # Deduplication of UMI
            output_prefix = "umi-deduplicated"
            print(f"UMI deduplication for {job_id} is now running.")
            subprocess.run((  # Run deduplication
                f"cd {job_dir}; "
                f"{shutil.which('umi_tools')} dedup "  # Define which umi_tools installation to use
                f"-I {after_sort} "  # Input
                f"--output-stats={output_prefix} "
                f"--paired "
                f"--unpaired-reads discard "
                f"-S {output_prefix}.bam "
                f"> report_{output_prefix}.log"
            ), shell=True)

            subprocess.run((  # Convert to sam file
                f"cd {job_dir}; "
                f"{shutil.which('/usr/bin/samtools')} view -h "  # Define which samtools installation to use
                f"-o {output_prefix}.sam {output_prefix}.bam"  # Output & Input
            ), shell=True)

            self.jobs[job_id]["process_genome_alignment"] = os.path.join(job_dir, f"{output_prefix}.sam")

        elif jobbing["sequencing_method"] in ["single", "paired_linking"]:

            print(f"Genome alignment for {job_id} is now running.")
            prefix = "genome_alignment_"
            subprocess.run((
                f"cd {job_dir}; " 
                f"{shutil.which('STAR')} "  # Define which star installation to use
                f"--runThreadN {self.cpu - 2} "  # Define how many core to be used. All cores are now using
                f"--genomeDir {self.index_directories['index_dna']} "
                f"--readFilesIn {input_fastq_fasta} "
                # All parameters were inherited from Mati-Kai's pipeline.
                "--outFilterMultimapNmax 1 "
                # "--peOverlapNbasesMin 6 " It is for paired end sequencing
                # "--peOverlapMMp 0.1 " It is for paired end sequencing
                "--outFilterType BySJout "
                "--alignIntronMin 5 "
                f"--outFileNamePrefix {os.path.join(job_dir, prefix)} "
                "--outReadsUnmapped Fastx "
                "--outSAMtype BAM SortedByCoordinate "
                "--outSAMattributes All XS "
                "--quantMode GeneCounts "
                "--twopassMode Basic "
                "> report_genome_alignment.log"
            ), shell=True)
            raw_bam = os.path.join(job_dir, prefix + "Aligned.sortedByCoord.out.bam")

            after_sort = os.path.splitext(raw_bam)[0] + '_sorted.bam'
            subprocess.run((  # Sort bam file and Index
                f"cd {job_dir}; "
                f"{shutil.which('/usr/bin/samtools')} sort {raw_bam} -o {after_sort}; "
                f"{shutil.which('/usr/bin/samtools')} index {after_sort}"
            ), shell=True)

            # Deduplication of UMI
            output_prefix = "umi-deduplicated"
            print(f"UMI deduplication for {job_id} is now running.")
            subprocess.run((  # Run deduplication
                f"cd {job_dir}; "
                f"{shutil.which('umi_tools')} dedup "  # Define which umi_tools installation to use
                f"-I {after_sort} "  # Input
                f"--output-stats={output_prefix} "
                f"-S {output_prefix}.bam "
                f"> {output_prefix}.log "
                f"> report_{output_prefix}.log"
            ), shell=True)

            subprocess.run((  # Convert to sam file
                f"cd {job_dir}; "
                f"{shutil.which('/usr/bin/samtools')} view -h "  # Define which samtools installation to use
                f"-o {output_prefix}.sam {output_prefix}.bam"  # Output & Input
            ), shell=True)

            self.jobs[job_id]["process_genome_alignment"] = os.path.join(job_dir, f"{output_prefix}.sam")

    def linking_pairs(self):
        for job_id in self.jobs:
            if 3 in self.jobs[job_id]["processes"]:
                self.linking_pairs_one_job(job_id)
            else:
                self.jobs[job_id]["process_linking_pairs"] = self.jobs[job_id]["process_cleanup"]

    def linking_pairs_one_job(self, job_id):

        jobbing = self.jobs[job_id]
        job_dir = jobbing["processes_dirs"][3]  # 1'i açıkla
        temp_fastq = self.jobs[job_id]["process_cleanup"]
        assert JobList.process_map[3].endswith("LinkingPairs")
        sam_path = os.path.join(job_dir, "prealignment.sam")
        print(f"Prealignment to link pairs for {job_id} is now running.")
        subprocess.run((
            f"cd {job_dir}; "  
            f"{shutil.which('bowtie2')} "  # Run Bowtie2 module
            "-D 40 -R 6 -N 0 -L 15 -i S,1,0.50 "  # Alignment effort and sensitivity. It is now very high.
            "-I20 -X120 "  # Search only those that has 20-120 nt. 
            "--score-min G,20,5.5 "  # Min score lowered
            "--ma 3 --mp 5,1 "  # ma bonus increased
            "--no-unal "  # To suppress the non-aligned reads
            f"-p{self.cpu - 2} "  # Number of core to use
            "--no-discordant "  # Filter pairs does not obey orientation/length constraints 
            "--no-mixed "  # Do not search for individual pairs if one in a pair does not align.
            "--local "  # Indicate the alignment is local. Soft clipping will be applied.
            "--time "  # Print the wall-clock time required to load the index files and align the reads. 
            f"-x {self.index_directories['index_cdna']} "  # Index directory with the base name
            "-q "  # Specifies the inputs are in fastq format
            f"-1 {temp_fastq[0]} "  # Read 1
            f"-2 {temp_fastq[1]} "  # Read 2
            f"-S {sam_path} "  # Output sam file
            "2> report_prealignment.log"
        ), shell=True)

        # Take the transcript sequences into random access memory
        fasta_transcriptome = self.org_db.get_db("cdna")
        transcript_seqs = dict()  # Create a dictionary with transcript id as keys
        with open(fasta_transcriptome, "r") as handle:  # Use previously filtered fasta file
            for record in SeqIO.parse(handle, "fasta"):
                transcript_seqs[record.id] = str(record.seq)  # Sequence as the values

        output_fasta = os.path.join(job_dir, "footprints.fasta")

        # Processing the SAM file
        print("Processing the SAM file.")
        with open(output_fasta, "w") as output_handle:  # Open output fasta file
            popup_dict = dict()
            with pysam.AlignmentFile(sam_path, "r") as sam_handle:  # Open sam file
                sam_iterator = sam_handle.fetch()  # Get the iterator
                for e in sam_iterator:  # Iterate over the entries

                    # Check if the entry is mapped, paired, and pair is mapped
                    if not e.is_unmapped and e.is_paired and e.is_proper_pair and not e.mate_is_unmapped:

                        # Determine start and end positions, taking into account soft-clipping from the ends.
                        qn = re.search(r"^[^_]*", e.query_name).group()  # Remove UMI info if exists
                        is_entry_complete = False
                        if not e.is_reverse and qn not in popup_dict:
                            popup_dict[qn] = [e.reference_start, None]
                        elif not e.is_reverse:  # and qn in popup_dict
                            popup_dict[qn][0] = e.reference_start
                            is_entry_complete = True
                        elif e.is_reverse and qn not in popup_dict:
                            popup_dict[qn] = [None, e.reference_end]
                        elif e.is_reverse:  # and qn in popup_dict
                            popup_dict[qn][1] = e.reference_end
                            is_entry_complete = True
                        else:
                            raise Exception("Unexpected entry!")

                        # Print the result if the dict is fully completed.
                        if is_entry_complete:
                            start_seq, end_seq = popup_dict.pop(qn)
                            # Different reference_names for pairs is impossible
                            fp = transcript_seqs[e.reference_name][start_seq: end_seq]
                            # Write down the identifier and sequence to fasta file
                            output_handle.write(f">{e.query_name}\n{fp}\n")

        self.jobs[job_id]["process_linking_pairs"] = output_fasta

    def cleanup(self):
        for job_id in self.jobs:
            if 2 in self.jobs[job_id]["processes"]:
                self.cleanup_one_job(job_id)
            else:
                self.jobs[job_id]["process_cleanup"] = self.jobs[job_id]["process_umitools"]

    def cleanup_one_job(self, job_id):

        jobbing = self.jobs[job_id]
        job_dir = jobbing["processes_dirs"][2]  # 2'i açıkla
        temp_paths = self.jobs[job_id]["process_umitools"]
        assert JobList.process_map[2].endswith("Cleanup")

        if jobbing["sequencing_method"] in ["paired", "paired_linking"]:
            output_path = os.path.join(job_dir, "Read%_norRNA.fastq")
            print(f"rRNA removal for {job_id} is now running.")
            subprocess.run((
                f"cd {job_dir}; "  # Change the directory to the index directory
                f"{shutil.which('bowtie2')} "  # Run Bowtie2 module
                f"-p{self.cpu - 2} "  # Number of core to use
                "--no-mixed "  # Do not search for individual pairs if one in a pair does not align.
                "-I20 -X120 "  # Default -I=0, -X=500. Since I will disregard X>120 and I<20 in the link-pairing module
                "--time "  # Print the wall-clock time required to load the index files and align the reads.
                "--score-min G,20,6 "  # Allow looser alignment. --ma 3 
                "--local --sensitive-local "  # Apply soft clipping when aligning. Default sensitivity.
                "-k1 "  # We are not interested in best alignment as long as it aligns somewhere in the indexed fasta.
                f"-x {self.index_directories['index_rrna']} "  # Index directory with the base name
                "-q "  # Indicates the inputs are fastq
                f"-1 {temp_paths[0]} "  # Read 1
                f"-2 {temp_paths[1]} "  # Read 2
                f"--un-conc {output_path} "  # Output fastq file, Contains all reads which did not aligned RNAs. 
                "-S /dev/null "  # Discard alignment sam file /dev/null
                f"2> report_cleanup.log"
            ), shell=True)

            self.jobs[job_id]["process_cleanup"] = [os.path.join(job_dir, j)
                                                    for j in ["Read1_norRNA.fastq", "Read2_norRNA.fastq"]]

        elif jobbing["sequencing_method"] in ["single"]:
            output_path = os.path.join(job_dir, "Read1_norRNA.fastq")
            print(f"rRNA removal for {job_id} is now running.")
            subprocess.run((
                f"cd {job_dir}; "  # Change the directory to the index directory
                f"{shutil.which('bowtie2')} "  # Run Bowtie2 module
                f"-p{self.cpu - 2} "  # Number of core to use
                "--time "  # Print the wall-clock time required to load the index files and align the reads.
                "--score-min G,20,6 "  # Allow looser alignment.  --ma 3 
                "--local --sensitive-local "  # Apply soft clipping when aligning. Default sensitivity.
                "-k1 "  # We are not interested in best alignment as long as it aligns somewhere in the indexed fasta.
                f"-x {self.index_directories['index_rrna']} "  # Index directory with the base name
                "-q "  # Indicates the inputs are fastq
                f"{temp_paths} "  # Read 1
                f"--un {output_path} "  # Output fastq file, Contains all reads which did not aligned RNAs. 
                "-S /dev/null "  # Discard alignment sam file /dev/null
                f"2> report_rnaremove.txt"
            ), shell=True)

            self.jobs[job_id]["process_cleanup"] = output_path

    def preprocessing(self):
        # explain1
        job_with_preprocessing = [job_id for job_id in self.jobs if 1 in self.jobs[job_id]["processes"]]
        for job_id in job_with_preprocessing:
            self.preprocessing_cutadapt(job_id)
        self.preprocessing_umitools_multiprocessing(job_with_preprocessing)
        for job_id in [job_id for job_id in self.jobs if 1 not in self.jobs[job_id]["processes"]]:
            self.jobs[job_id]["process_cutadapt"] = self.jobs[job_id]["input_fastq"]
            self.jobs[job_id]["process_umitools"] = self.jobs[job_id]["input_fastq"]

    def preprocessing_umitools_multiprocessing(self, job_id_list):
        assert len(self.jobs) <= self.cpu - 2
        executor = multiprocessing.Pool(len(self.jobs))
        result = executor.map(self.preprocessing_umitools, job_id_list)
        executor.terminate()
        executor.join()
        for job_id, final_paths in zip(job_id_list, result):
            self.jobs[job_id]["process_umitools"] = final_paths

    def preprocessing_umitools(self, job_id):

        jobbing = self.jobs[job_id]
        job_dir = jobbing["processes_dirs"][1]  # 1'i açıkla
        temp_paths = self.jobs[job_id]["process_cutadapt"]
        assert JobList.process_map[1].endswith("Preprocessing")

        if jobbing["sequencing_method"] in ["paired", "paired_linking"]:
            final_paths = [f"{i}_no-adapt_umi-aware.fastq.gz" for i in ["read1", "read2"]]
            final_paths = [os.path.join(job_dir, i) for i in final_paths]
            pattern1, pattern2 = jobbing["pattern_umi1"], jobbing["pattern_umi2"]
            print(f"Umitool for {job_id} is now running.")
            subprocess.run((
                f"cd {job_dir}; "  # Change the directory to the index directory
                f"{shutil.which('umi_tools')} extract "  # Define which extract installation to use
                "--extract-method=regex "
                f"--bc-pattern='{pattern1}' "  # Barcode pattern
                f"--bc-pattern2='{pattern2}' "  # Barcode pattern for paired reads"
                f"-I {temp_paths[0]} -S {final_paths[0]} "  # Input and output for read 1
                f"--read2-in={temp_paths[1]} --read2-out={final_paths[1]} "  # Input and output for read 2
                f"--log=umi_tools.log "  # Log the results
                f"> 'report_umitool.log'"
            ), shell=True)

            return final_paths

        elif jobbing["sequencing_method"] in ["single"]:
            final_path = os.path.join(job_dir, "read_1_no-adapt_umi-aware.fastq.gz")
            umitools_pattern = jobbing["pattern_umi"]
            print(f"Umitool for {job_id} is now running.")
            subprocess.run((
                f"cd {job_dir}; "  # Change the directory to the index directory
                f"{shutil.which('umi_tools')} extract "  # Define which extract installation to use
                "--extract-method=regex "
                f"--bc-pattern='{umitools_pattern}' "  # Barcode pattern
                f"-I {temp_paths} -S {final_path} "  # Input and output for read 1
                f"--log=umi_tools.log "  # Log the results
                f"> 'report_umitool.log'"
            ), shell=True)

            return final_path

    def preprocessing_cutadapt(self, job_id):

        jobbing = self.jobs[job_id]
        job_dir = jobbing["processes_dirs"][1]  # 1'i açıkla
        input_fastq = self.jobs[job_id]["input_fastq"]
        assert JobList.process_map[1].endswith("Preprocessing")

        if jobbing["sequencing_method"] in ["paired", "paired_linking"]:
            # Outputs for the first run
            temp_paths = [f"read1_cutadapt_temp.fastq.gz", f"read2_cutadapt_temp.fastq.gz"]
            temp_paths = [os.path.join(job_dir, i) for i in temp_paths]
            # Create flag if adapter is provided
            read1_adapter = f"-a {jobbing['adapter1']}" if jobbing['adapter1'] else ""
            read2_adapter = f"-A {jobbing['adapter2']}" if jobbing['adapter2'] else ""
            print(f"Cutadapt for {job_id} is now running.")
            subprocess.run((
                    f"cd {job_dir}; "  # Change the directory to the index directory
                    f"{shutil.which('cutadapt')} "  # Define which cutadapt installation to use
                    f"--cores={self.cpu - 2} "  # Define how many core to be used. All cores are now using
                    "--match-read-wildcards "
                    f"-q20 -m23 "  # Settings for quality and minimum length to cut the adapter sequence
                    # todo: explain why: check the photo taken in 12 february
                    # "--discard-untrimmed "
                    # "-O6 "  # Require minimum length overlap between read and adapter for an adapter to be found.
                    f"{read1_adapter} {read2_adapter}".strip() + " "  # Adapter flags. No flanking white space allowed
                    f"-o {temp_paths[0]} -p {temp_paths[1]} "  # Path to output trimmed sequences
                    f"{input_fastq[0]} {input_fastq[1]} "  # Input file
                    f"1> 'report_cutadapt_temp.log'"
            ), shell=True)

            self.jobs[job_id]["process_cutadapt"] = temp_paths

        elif jobbing["sequencing_method"] in ["single"]:
            temp_path = os.path.join(job_dir, "read_1_cutadapt_temp.fastq.gz")  # Outputs for the first run
            print(f"Cutadapt for {job_id} is now running.")
            subprocess.run((
                f"cd {job_dir}; "  # Change the directory to the index directory
                f"{shutil.which('cutadapt')} "  # Define which cutadapt installation to use
                f"--cores={self.cpu - 2} "  # Define how many core to be used. All cores are now using
                "--match-read-wildcards "
                f"-q20 -m23 "  # Settings for quality and minimum length to cut the adapter sequence
                # If we don't below two, UMI tool definitely malfunction
                "--discard-untrimmed "
                "-O6 "  # Require minimum length overlap between read and adapter for an adapter to be found.
                f"-a {jobbing['adapter']} "  # The adapter flags. No flanking white space allowed
                f"-o {temp_path} "  # Path to output trimmed sequences
                f"{input_fastq} "  # Input file
                f"1> 'report_cutadapt_temp.log'"
            ), shell=True)

            self.jobs[job_id]["process_cutadapt"] = temp_path

    def create_output_tree(self):
        dir_run_identifier = create_dir(self.data_repo_dir, self.run_identifier)
        for job_id in self.jobs:
            dir_job = create_dir(dir_run_identifier, job_id)
            self.jobs[job_id]["processes_dirs"] = {process: create_dir(dir_job, JobList.process_map[process])
                                                   for process in self.jobs[job_id]["processes"]}

    def create_index(self):
        dir_base = create_dir(self.temp_repo_dir, self.organism)
        dir_base_index = create_dir(dir_base, "index")
        dir_cdna = create_dir(dir_base_index, "cdna")
        dir_dna = create_dir(dir_base_index, "dna")
        dir_rrna = create_dir(dir_base_index, "rrna")
        index_name_cdna = f"{self.organism}_cdna"
        index_name_rrna = f"{self.organism}_rrna"

        if not self.create_index_search_metadata(dir_cdna):
            print("Indexing for cDNA is being calculated.")
            temp_cdna_fasta = self.org_db.get_db("cdna")
            subprocess.run((
                f"cd {dir_cdna}; "  # Change the directory to the index directory
                f"{shutil.which('bowtie2-build')} "  # Name of the function
                f"--threads {self.cpu - 2} "  # Number of threads to be used.
                f"{temp_cdna_fasta} "  # Input file. -f is to indicate the file is in fasta format
                f"{index_name_cdna} "  # The basename of the index files to write
                f"> report_{self.organism}_cdna.log"
            ), shell=True)
            metadata_index = get_files_metadata(dir_cdna)  # Write the file info
            joblib.dump(metadata_index, os.path.join(dir_cdna, ".metadata.joblib"))

        if not self.create_index_search_metadata(dir_rrna):
            print("Indexing for rRNA is being calculated.")
            temp_rrna_fasta = self.org_db.get_db("rrna")
            subprocess.run((
                f"cd {dir_rrna}; "  # Change the directory to the index directory
                f"{shutil.which('bowtie2-build')} "  # Name of the function
                f"--threads {self.cpu - 2} "  # Number of threads to be used.
                f"{temp_rrna_fasta} "  # Input file. -f is to indicate the file is in fasta format
                f"{index_name_rrna} "  # The basename of the index files to write
                f"> report_{self.organism}_rrna.log"
            ), shell=True)
            metadata_index = get_files_metadata(dir_rrna)  # Write the file info
            joblib.dump(metadata_index, os.path.join(dir_rrna, ".metadata.joblib"))

        if not self.create_index_search_metadata(dir_dna):
            # Manual says: "In most cases, the default value of 100 will work as well as the ideal value."
            read_length_minus_1 = 100
            print("Indexing for DNA is being calculated.")
            temp_dna_fasta = self.org_db.get_db("dna")
            temp_gtf_fasta = self.org_db.get_db("gtf")
            subprocess.run((
                f"cd {dir_dna}; "  # Change the directory to the index directory
                f"{shutil.which('STAR')} "  # Define which star installation to use
                f"--runThreadN {self.cpu - 2} "  # Define how many core to be used. All cores are now using
                "--runMode genomeGenerate "
                f"--genomeDir {dir_dna} "  # Directory to save the files
                f"--genomeFastaFiles {temp_dna_fasta} "  # Specifies FASTA file with the genome reference sequences
                f"--sjdbGTFfile {temp_gtf_fasta} "  # Specifies the file with annotated transcripts in the GTF format
                f"--sjdbOverhang {read_length_minus_1} "  # The length of the sequence around the annotated junction
                f"> report_{self.organism}_dna.log"
            ), shell=True)
            metadata_index = get_files_metadata(dir_dna)  # Write the file info
            joblib.dump(metadata_index, os.path.join(dir_dna, ".metadata.joblib"))

        return {"index_cdna": os.path.join(dir_cdna, index_name_cdna),
                "index_rrna": os.path.join(dir_rrna, index_name_rrna),
                "index_dna": dir_dna}

    @staticmethod
    def create_index_search_metadata(dir_path):
        # Check if the index file
        try:
            index_metadata_path = os.path.join(dir_path, ".metadata.joblib")
            assert os.path.isfile(index_metadata_path) and os.access(index_metadata_path, os.R_OK)
            metadata_index_previously = joblib.load(index_metadata_path)
            metadata_index = get_files_metadata(dir_path)
            assert metadata_index_previously == metadata_index
            return True
        except AssertionError:
            return False


class OrganismDatabase:

    ID_MAP = {"homo_sapiens": 9606, "mus_musculus": 10090}  # todo

    def __init__(self, organism, ensembl_release, temp_repo_dir):
        assert organism in ["homo_sapiens", "mus_musculus"]  # todo
        # todo: also add "homo_sapiens_refseq"

        self.ensembl_release = ensembl_release
        self.organism = organism
        self.organism_id = OrganismDatabase.ID_MAP[self.organism]
        self.temp_repo_dir = temp_repo_dir
        self.repository = create_dir(self.temp_repo_dir, self.organism)

        base_temp = f"ftp://ftp.ensembl.org/pub/release-{ensembl_release}"
        
        if self.organism == "homo_sapiens":
            # Genome GTF
            gtf_temp = f"gtf/homo_sapiens/Homo_sapiens.GRCh38.{self.ensembl_release}.chr_patch_hapl_scaff.gtf.gz"
            self.gtf = os.path.join(base_temp, gtf_temp)
            # Genome GFF3
            gff3_temp = f"gff3/homo_sapiens/Homo_sapiens.GRCh38.{self.ensembl_release}.chr_patch_hapl_scaff.gff3.gz"
            self.gff3 = os.path.join(base_temp, gff3_temp)
            # Genome DNA fasta
            dna_temp = "fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
            self.dna = os.path.join(base_temp, dna_temp)
            # Transcriptome DNA fasta
            cdna_temp = "fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
            self.cdna = os.path.join(base_temp, cdna_temp)
        
        elif self.organism == "mus_musculus":
            # Genome GTF
            gtf_temp = f"gtf/mus_musculus/Mus_musculus.GRCh38.{self.ensembl_release}.chr_patch_hapl_scaff.gtf.gz"
            self.gtf = os.path.join(base_temp, gtf_temp)
            # Genome GFF3
            gff3_temp = f"gff3/mus_musculus/Mus_musculus.GRCh38.{self.ensembl_release}.chr_patch_hapl_scaff.gff3.gz"
            self.gff3 = os.path.join(base_temp, gff3_temp)
            # Genome DNA fasta
            dna_temp = "fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
            self.dna = os.path.join(base_temp, dna_temp)
            # Transcriptome DNA fasta
            cdna_temp = "fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"
            self.cdna = os.path.join(base_temp, cdna_temp)

        rrna_base_temp = "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release"
        self.rrna_raw_fasta = os.path.join(rrna_base_temp, "sequences/by-database/ena.fasta")
        self.rrna_raw_information = os.path.join(rrna_base_temp, "id_mapping/database_mappings/ena.tsv")

    def get_db(self, db):
        if db == "rrna":
            return self._filter_rrna()
        db_url = eval(f"self.{db}")
        output_path_compressed = os.path.join(self.repository, os.path.basename(db_url))
        output_path_uncompressed = os.path.splitext(output_path_compressed)[0]
        if not os.access(output_path_uncompressed, os.R_OK) or not os.path.isfile(output_path_uncompressed):
            print(f"Downloading from the server for {db}:\n{db_url}")
            if not os.access(output_path_compressed, os.R_OK) or not os.path.isfile(output_path_compressed):
                subprocess.run(f"cd {self.repository}; curl -L -O --silent {db_url}", shell=True)
            subprocess.run(f"cd {self.repository}; gzip -d -q {output_path_compressed}", shell=True)
        return output_path_uncompressed

    def _download_rna_central(self, db_url):
        output_path = os.path.join(self.repository, os.path.basename(db_url))
        if not os.access(output_path, os.R_OK) or not os.path.isfile(output_path):
            print(f"Downloading from the server for rrna:\n{db_url}")
            subprocess.run(f"cd {self.repository}; curl -L -O --silent {db_url}", shell=True)
        return output_path

    def _filter_rrna(self):
        db_info = self._download_rna_central(self.rrna_raw_information)
        db_fasta = self._download_rna_central(self.rrna_raw_fasta)
        db_fasta_filtered = os.path.splitext(db_fasta)[0] + "_only_rrna.fasta"
        if not os.access(db_fasta_filtered, os.R_OK) or not os.path.isfile(db_fasta_filtered):
            print(f"Filtering rRNAs for {self.organism}")
            rrna_ids = set()
            with open(db_info, "r") as handle_info:
                for line in handle_info:
                    line = line.split()
                    if line[-1] == "rRNA":
                        rrna_ids.add(line[0])
            str_id = str(self.organism_id)
            with open(db_fasta, "rt") as handle_fasta, open(db_fasta_filtered, "wt") as handle_output:
                for record in SeqIO.parse(handle_fasta, "fasta"):
                    split_record_id = record.id.split("_")
                    if str_id == split_record_id[1] and split_record_id[0] in rrna_ids:
                        SeqIO.write(record, handle_output, "fasta")
        return db_fasta_filtered


class Col:
    H = '\033[95m\033[1m'  # Header
    B = '\033[94m'  # Blue
    C = '\033[96m'  # Cyan
    G = '\033[92m'  # Green
    W = '\033[93m'  # Warning
    F = '\033[91m'  # Fail
    E = '\033[0m'  # End
    D = '\033[1m'  # Bold
    U = '\033[4m'  # Underline


def create_dir(*args):
    """
    Creates a directory if not exist.
    :param args: Strings to use in os.path.join
    :return: Path of created file
    """
    dir_path = os.path.join(*args)
    if not os.access(dir_path, os.W_OK) or not os.path.isdir(dir_path):  # Create directory if not exist
        os.mkdir(dir_path)
    return dir_path


def check_exist_package(packages):
    for pkg in packages:
        if isinstance(pkg, str) and not shutil.which(pkg):
            raise ModuleNotFoundError(pkg)


def scan_tree(mother_directory):
    """Recursively yield DirEntry objects for given directory."""
    for entry in os.scandir(mother_directory):
        if entry.name.startswith('.'):
            continue
        if entry.is_dir(follow_symlinks=False):
            yield from scan_tree(entry.path)  # see below for Python 2.x
        else:
            yield entry


def get_files_metadata(mother_directory):
    scan_tree_iterator = scan_tree(mother_directory)
    output = list()
    for i in scan_tree_iterator:
        if not i.name.startswith('.') and i.is_file() and not i.is_symlink():
            stats = i.stat()
            output.append([i.path.split(mother_directory)[1], stats.st_size, stats.st_mtime, stats.st_ctime])
    return sorted(output)


def check_directory(paths):
    for path in paths:
        if os.path.isdir(path) and os.access(path, os.W_OK) and os.access(path, os.R_OK):
            return None
        else:
            raise NotADirectoryError(path)


def check_file(file_path, extensions):
    return os.path.isfile(file_path) and \
           os.access(file_path, os.R_OK) and \
           any([file_path.endswith(i) for i in extensions])


if __name__ == '__main__':
    main()


# End
