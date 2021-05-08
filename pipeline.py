#!/usr/bin/env python

import argparse
import datetime
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import math
import psutil

import joblib
import pysam
from Bio import SeqIO


# Authorship information
__author__ = "Kemal İnecik"
__license__ = "GPLv3"
__version__ = "0.2"
__maintainer__ = "Kemal İnecik"
__email__ = "k.inecik@gmail.com"
__status__ = "Development"


def main():
    """
    If running the module (the source file) as the main program, this function will be called.
    :return: None
    """
    a = argument_parser()  # Get the command line arguments
    # Create the job dictionary by using the introductory txt file.
    job_list = JobList(a.filepath)
    # Make sure the process will be run with correct settings by making it confirmed by the user.
    job_list.confirm_job_list()
    print(f"{Col.H}Operation started.{os.linesep}{Col.E}")  # Print that the processing has just started.
    # Create a controller object, which is the main functional object processing the riboseq data.
    controller = Controller(a.identifier, a.organism, a.ensembl_release, a.cpu, a.temp,
                            a.output, a.assignment, job_list.jobs)
    controller.start_processing()  # Start processing by explicitly calling it.
    print(f"{Col.H}Operation successfully ended.{os.linesep}{Col.E}")  # Print that the processing has just ended.


def argument_parser() -> argparse.Namespace:
    """
    Extract the information given by the user via command line. It also contains the explanation of the flags and
    the program itself
    :return: NameSpace object containing parsed arguments.
    """

    parser = argparse.ArgumentParser(prog="pipeline.py",  # Name of the program
                                     description="The app description.")  # todo

    # Introduce arguments, they are already self explanatory.
    parser.add_argument("-r", type=str, dest="identifier",
                        help="the identifier for this run, which will be the directory name for all outputs under "
                             "provided main output directory. If skipped, the default is a string containing "
                             "date and time of the program start.",
                        required=False, default=datetime.datetime.now().strftime("Run_on_%Y_%m_%d_at_%H_%M_%S"))

    parser.add_argument("-a", type=str, dest="organism",
                        help="organism of interest for the analysis. If skipped, the default value is homo_sapiens.",
                        required=False, default="homo_sapiens",
                        choices=["homo_sapiens", "mus_musculus", "saccharomyces_cerevisiae", "escherichia_coli"],
                        )

    parser.add_argument("-e", type=int, dest="ensembl_release",
                        help="ensembl version to be used. Ignored if the 'organism' does not necessitates Ensembl "
                             "release information. Default value is 102. For E. coli, it must be 48.",
                        required=False, default=102)

    parser.add_argument("-c", type=int, dest="cpu",
                        help="number of cpu cores to be used. Default value is maximum minus eight.",
                        required=False, default=multiprocessing.cpu_count() - 8)

    parser.add_argument("-s", type=int, dest="assignment",
                        help="select 3' or 5' assignment for Julia script.",
                        required=False, default=3,
                        choices=[3, 5])

    parser.add_argument("-f", type=str, required=True, dest="filepath",
                        help="path of the txt file which contains the task list.")

    parser.add_argument("-t", type=str, required=True, dest="temp",
                        help="absolute path of the directory to be used for temporary files such as genome indexes.")

    parser.add_argument("-o", type=str, required=True, dest="output",
                        help="absolute path of the directory to be used for output files.")

    return parser.parse_args()  # Return the parsed result.


class JobList:
    """
    Processes the introductory txt file, raises errors if the input is not fine, creates a dictionary with the content,
    also makes the result confirmed by the user.
    """

    # Steps of the program. Changing the numbers here potentially makes whole pipeline to malfunction. If you will add
    # any additional step, please take this into consideration. Either add the new step after 5, or just edit whole
    # script properly.
    process_map = {1: "01_Preprocessing",
                   2: "02_Cleanup",
                   3: "03_LinkingPairs",
                   4: "04_GenomeAlignment",
                   5: "05_JuliaAssignment"}
    # Default values are based on the data published on Matilde et. al 2021.
    # Adapter sequences
    default_adapter_single = "ATCGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    default_adapter1_paired = "ATCGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    default_adapter2_paired = None  # None, if no adapter trimming is required.
    # UMI patterns, the below patterns are regex strings.
    default_pattern_single = "^(?P<umi_1>.{2}).*(?P<umi_2>.{5})$"
    default_pattern1_paired = "^(?P<umi_1>.{2}).*"
    default_pattern2_paired = "^(?P<discard_2>.{5})(?P<umi_2>.{5}).*"
    # Default processing steps, the numbers corresponds to the steps of process_map dictionary above.
    default_processes_single = [1, 2, 4, 5]
    default_processes_paired_linking = [1, 2, 3, 4, 5]
    default_processes_paired = [1, 2, 4]
    # For now, the only valid input for raw processing files are fastq.gz. Other file formats were not tested.
    proper_input_extensions = ["fastq.gz"]

    def __init__(self, user_task_file: str):
        """
        Initiates the object.
        :param user_task_file: Absolute or relative path of introductory txt file
        """
        self.user_task_file = user_task_file
        self.jobs = self.read_job_list(self.user_task_file)

    @staticmethod
    def read_job_list(file_path: str) -> dict:
        """
        Processes the introductory txt file, raises errors if the input file is not in the allowed format, creates a
        dictionary with the file content.
        :param file_path: Absolute or relative path of introductory txt file
        :return: Job IDs as the keys relevant information as the values.
        """
        result = dict()  # Create a dictionary to populate with jobs
        with open(file_path, "r") as handle:  # Open the file in read only mode.
            # Read all the data in the file. Split the content first with '>>>' to separate out the jobs. Skip comments
            content = f"{os.linesep}".join([line.strip() for line in handle if not line.startswith('#')])

            for entry in [r.strip() for r in content.split(">>>") if r.strip()]:  # For each job
                # Separate out the lines in a given job by new line character.
                # If it starts with '#' just ignore this line
                line = [r.strip() for r in entry.split(os.linesep) if r.strip()]
                print_error_entry = Col.B + os.linesep.join(line) + Col.E  # A string to be printed in case of an error.
                # Do not allow whitespace inside an information. Note that lines were already striped above.
                if not all([len(i.split()) == 1 for i in line]):
                    raise ValueError(f"No whitspace is allowed:{os.linesep}{print_error_entry}")
                # Make sure the second line (first starts with '>>>') is one of below list.
                if line[1] not in ["single", "paired_linking", "paired"]:
                    raise ValueError(f"Improper sequencing method:{os.linesep}{print_error_entry}")
                # Make sure input raw file is exists, readable, and with one of the allowed extensions.
                if not check_file(line[2], JobList.proper_input_extensions):
                    raise ValueError(f"Improper file or file extension method:{os.linesep}{print_error_entry}")
                # Make sure that input raw file path is absolute.
                if not os.path.isabs(line[2]):
                    raise ValueError(f"File paths have to be absolute:{os.linesep}{print_error_entry}")
                # Make sure input raw file is exists, readable, and with one of the allowed extensions.
                if line[1].startswith("paired") and not check_file(line[3], JobList.proper_input_extensions):
                    raise ValueError(f"Improper file or file extension method:{os.linesep}{print_error_entry}")
                # Make sure that input raw file path is absolute.
                if line[1].startswith("paired") and not os.path.isabs(line[3]):
                    raise ValueError(f"File paths have to be absolute:{os.linesep}{print_error_entry}")

                # If the job is for single end sequencing. Initiate the job dictionary with default elements.
                if line[1] == "single":
                    temp_dict = {"sequencing_method": line[1],
                                 "input_fastq": line[2],
                                 "adapter": JobList.default_adapter_single,
                                 "pattern_umi": JobList.default_pattern_single,
                                 "processes": JobList.default_processes_single}
                else:  # If the job is for paired end sequencing. Initiate the job dictionary with default elements.
                    # The 'paired_linking' and 'paired' needs to be different in processes parameter, since
                    # paired_linking requires an extra LinkingPairs step.
                    process_temp = JobList.default_processes_paired if line[1] == "paired" \
                        else JobList.default_processes_paired_linking
                    temp_dict = {"sequencing_method": line[1],
                                 "input_fastq": [line[2], line[3]],
                                 "adapter1": JobList.default_adapter1_paired,
                                 "adapter2": JobList.default_adapter2_paired,
                                 "pattern_umi1": JobList.default_pattern1_paired,
                                 "pattern_umi2": JobList.default_pattern2_paired,
                                 "processes": process_temp}
                # If the user want to change adapters, umi, or processes, the lines after 3 for single end sequencing,
                # and the lines after 4 for paired end sequencing will be populated.
                non_default_entries = line[3:] if line[1] == "single" else line[4:]
                # Allowed keys for the non default entries
                if line[1] == "single":
                    possible_settings = ["adapter", "pattern_umi", "processes"]
                else:
                    possible_settings = ["adapter1", "adapter2", "pattern_umi1", "pattern_umi2", "processes"]
                for nde in non_default_entries:  # For each non default entry
                    # First split with '=', and strip the output.
                    des = [r.strip() for r in nde.split("=") if r.strip()]
                    if len(des) != 2 or des[0] not in possible_settings:  # Ignore if setting is not with allowed keys.
                        # If there is two '=' there will be more than 2 elements after split. Alternatively, if there
                        # is no '=' or one side of the '=' is not populated, there will be one element.
                        print(f"{Col.W}Manual setting ignored for {line[0]}:{os.linesep}{des}{Col.E}")
                    elif des[0] == "processes":  # If the user wants to customize processes.
                        # Split by ',' and convert everything into integers. Error will be raised if it is not possible.
                        temp_dict[des[0]] = sorted([int(r.strip()) for r in des[1].split(',') if r.strip()])
                    else:  # For other customizations.
                        # Just get the value directly if the value is not 'None'. If 'None', set the value to None.
                        temp_dict[des[0]] = des[1] if des[1] != "None" else None
                # Make sure that all entries are unique
                assert line[0] not in result, f"Duplicated key for {line[0]}"
                # Add the dictionary for the individual job to the dictionary of all jobs.
                result[line[0]] = temp_dict
        return result  # Return the final dictionary.

    def confirm_job_list(self):
        """
        Method to make sure the processing the txt file is correct. It prints out the jobs and request confirmation
        by the user. If the user types anything except 'enter', the program will be aborted.
        :return: None
        """
        # Print the instruction.
        print(f"{Col.H}Please confirm the jobs.{os.linesep}"
              f"Press enter to confirm, type anything to stop.{os.linesep}"
              f"There are {len(self.jobs)} jobs:{Col.E}{os.linesep}")

        for job_id in self.jobs:  # For all jobs in the job list.
            the_job = self.jobs[job_id]  # Get the job information dictionary.
            # Convert the processing into str with comma as the delimiter to be printed out
            processes_pre = ", ".join([str(i) for i in the_job['processes']])

            if the_job["sequencing_method"] == "single":  # If the job is a single end sequencing job.
                # If the default values are used, just add 'Default'. Otherwise, print the information as it is.
                adapter_temp = the_job["adapter"] if the_job["adapter"] != JobList.default_adapter_single \
                    else f"Default (\"{the_job['adapter']}\")"
                pattern_temp = the_job["pattern_umi"] if the_job["pattern_umi"] != JobList.default_pattern_single \
                    else f"Default (\"{the_job['pattern_umi']}\")"
                processes_temp = processes_pre if the_job["processes"] != JobList.default_processes_single \
                    else f"Default (\"{processes_pre}\")"
                # Prepare a string to be printed out in the command line
                print_string = (
                    f"{Col.H}> {job_id}{Col.C}{os.linesep}"  # Job ID
                    f"Sequencing Method : {the_job['sequencing_method']}{os.linesep}" 
                    f"Read              : {the_job['input_fastq']}{os.linesep}"
                    f"Adapter           : {adapter_temp}{os.linesep}"
                    f"Pattern UMI       : {pattern_temp}{os.linesep}"
                    f"Processes         : {processes_temp}{Col.E}{os.linesep}"
                )

            else:  # If the job is a paired end sequencing job.
                # If the default values are used, just add 'Default'. Otherwise, print the information as it is.
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
                # Prepare a string to be printed out in the command line
                print_string = (
                    f"{Col.H}> {job_id}{Col.C}{os.linesep}"  # Job ID
                    f"Sequencing Method : {the_job['sequencing_method']}{os.linesep}"
                    f"Read 1            : {the_job['input_fastq'][0]}{os.linesep}"
                    f"Read 2            : {the_job['input_fastq'][1]}{os.linesep}"
                    f"Adapter 1         : {adapter1_temp}{os.linesep}"
                    f"Adapter 2         : {adapter2_temp}{os.linesep}"
                    f"Pattern UMI 1     : {pattern1_temp}{os.linesep}"
                    f"Pattern UMI 2     : {pattern2_temp}{os.linesep}"
                    f"Processes         : {processes_temp}{Col.E}{os.linesep}"
                )

            confirm = input(print_string)  # Request the user's response
            if confirm != "":  # If not pressed anything except directly 'enter'
                print(f"{Col.F}Process terminated.{Col.E}")
                sys.exit(1)  # Terminate the process with an error code 1.
            else:  # Otherwise just move to the next job.
                print(f"{Col.W}Confirmed.{os.linesep}{Col.E}")


class Controller:

    def __init__(self, run_identifier: str, organism: str, ensembl_release: int, cpu_cores: int,
                 temp_repo_dir: str, data_repo_dir: str, assign_from: int, jobs: dict):
        # First check if the directories for temp and output, make sure they exist and writable.
        check_directory([temp_repo_dir, data_repo_dir])
        # Make sure the below third party packages are installed.
        # Systems samtools will be used because the conda distribution raises some errors.
        check_exist_package(["cutadapt", "umi_tools", "bowtie2-build", "STAR", "bowtie2", "/usr/bin/samtools"])
        # Assign the parameters as object variables
        self.temp_repo_dir = temp_repo_dir
        self.data_repo_dir = data_repo_dir
        self.organism = organism
        self.ensembl_release = ensembl_release
        self.jobs = jobs
        self.run_identifier = run_identifier
        self.cpu = cpu_cores
        self.assign_from = assign_from
        # Make sure julia assignment script is in the same directory with this script.
        self.julia_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "julia_assignment.jl")
        assert check_file(self.julia_path, [".jl"]), f"JuliaAssignment script could not be found.{os.linesep}:" \
                                                     f"{self.julia_path}"
        # Create an OrganismDatabase object with designated parameters.
        self.org_db = OrganismDatabase(self.organism, self.ensembl_release, self.temp_repo_dir)
        # Call the method to create the folders under data_repo_dir to save the final and intermediate results.
        self.create_output_tree()
        # Call the method to create index files for genome, rRNA and transcriptome alignments. Note that only
        # pair_linking uses transcriptome indexes but this will be created anyway.
        self.index_directories = self.create_index()

    def is_already_calculated(self):
        pass  # todo

    def start_processing(self):
        """
        Runs the processing for all steps.
        :return: None
        """
        self.preprocessing()
        self.cleanup()
        self.linking_pairs()
        self.genome_alignment()
        self.julia_assignment()
        # Also save the Controller object as a joblib file, to be able to see running parameters in the future.
        joblib.dump(self, os.path.join(self.data_repo_dir, self.run_identifier, "pipeline_controller.joblib"))

    def julia_assignment(self):
        """
        For the jobs which is set to be processed, run julia assignment.
        :return: None
        """
        gff_path = self.julia_assignment_gff3_correct()  # Get GFF file path
        for job_id in self.jobs:  # Iterate over the jobs
            if 5 in self.jobs[job_id]["processes"]:  # If it is set to be processed by this step.
                self.julia_assignment_one_job(job_id, gff_path)  # Run the processing for this job.
        # The last step of the pipeline, no output path is returned.

    def julia_assignment_one_job(self, job_id, gff_path):
        """

        :param job_id:
        :param gff_path:
        :return:
        """
        jobbing = self.jobs[job_id]
        job_dir = jobbing["processes_dirs"][5]  # 5'i açıkla
        input_sam = self.jobs[job_id]["process_genome_alignment"]
        assert JobList.process_map[5].endswith("JuliaAssignment")

        run_and_check((
            f"cd {job_dir}; "  # Change the directory to the directory
            f"{shutil.which('julia')} {self.julia_path} "  # Which Julia installation to use and the script
            f"-g {gff_path} "  # Gff3 file. Removed of duplicated gene names
            f"-a {self.assign_from} "  # Assignment from 3'
            # "-u "  # Inherited from Matilde et al 2021. Removed here because umi-tool deduplication was already done.
            f"-o {job_dir} {input_sam}"  # Output file & Input file
        ), shell=True)  # The specified command will be executed through the shell.
        # The last step of the pipeline, no output path is returned.

    def julia_assignment_gff3_correct(self):
        gff_path = self.org_db.get_db("gff3")
        output_path = os.path.splitext(gff_path)[0] + "_renamed_duplicate_gene_names.gff3"

        try:

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
        except KeyError:
            return gff_path

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
        prefix = "genome_alignment_"

        if jobbing["sequencing_method"] in ["paired", "paired_linking"] \
                and (jobbing["pattern_umi1"] or jobbing["pattern_umi1"]):
            self.jobs[job_id]["is_umi_extracted"] = True
        elif jobbing["sequencing_method"] == "single" and jobbing["pattern_umi"]:
            self.jobs[job_id]["is_umi_extracted"] = True
        else:
            self.jobs[job_id]["is_umi_extracted"] = False

        available_memory = int(psutil.virtual_memory().available * 9 / 10)

        if jobbing["sequencing_method"] in ["paired"]:  # note it down, not paired-linking

            print(f"Genome alignment for {job_id} is now running.")
            run_and_check((
                f"cd {job_dir}; "  # Change the directory to the index directory
                f"{shutil.which('STAR')} "  # Define which star installation to use
                f"--runThreadN {self.cpu} "  # Define how many core to be used. All cores are now using
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
                f"--limitBAMsortRAM {available_memory} "  # If not set, error can be raised for some cases.
                "--outSAMattributes All XS "
                "--quantMode GeneCounts "
                "--twopassMode Basic "
                "> report_genome_alignment.log"
            ), shell=True)  # The specified command will be executed through the shell.
            raw_bam = os.path.join(job_dir, prefix + "Aligned.sortedByCoord.out.bam")

            if not jobbing["is_umi_extracted"]:
                raw_sam = os.path.join(job_dir, prefix + "Aligned.sortedByCoord.out.sam")
                run_and_check(f"cd {job_dir}; {shutil.which('/usr/bin/samtools')} view -h -o {raw_sam} {raw_bam}",
                              shell=True)  # The specified command will be executed through the shell.
                self.jobs[job_id]["process_genome_alignment"] = raw_sam
            else:
                after_sort = os.path.splitext(raw_bam)[0] + '_sorted.bam'
                run_and_check((  # Sort bam file and Index
                    f"cd {job_dir}; "
                    f"{shutil.which('/usr/bin/samtools')} sort {raw_bam} -o {after_sort}; "
                    f"{shutil.which('/usr/bin/samtools')} index {after_sort}"
                ), shell=True)  # The specified command will be executed through the shell.

                # Deduplication of UMI
                output_prefix = "umi-deduplicated"
                print(f"UMI deduplication for {job_id} is now running.")
                run_and_check((  # Run deduplication
                    f"cd {job_dir}; "
                    f"{shutil.which('umi_tools')} dedup "  # Define which umi_tools installation to use
                    f"-I {after_sort} "  # Input
                    f"--output-stats={output_prefix} "
                    f"--paired "
                    f"--unpaired-reads discard "
                    f"-S {output_prefix}.bam "
                    f"> report_{output_prefix}.log"
                ), shell=True)  # The specified command will be executed through the shell.

                run_and_check((  # Convert to sam file
                    f"cd {job_dir}; "
                    f"{shutil.which('/usr/bin/samtools')} view -h "  # Define which samtools installation to use
                    f"-o {output_prefix}.sam {output_prefix}.bam"  # Output & Input
                ), shell=True)  # The specified command will be executed through the shell.
                self.jobs[job_id]["process_genome_alignment"] = os.path.join(job_dir, f"{output_prefix}.sam")

        elif jobbing["sequencing_method"] in ["single", "paired_linking"]:

            print(f"Genome alignment for {job_id} is now running.")
            run_and_check((
                f"cd {job_dir}; " 
                f"{shutil.which('STAR')} "  # Define which star installation to use
                f"--runThreadN {self.cpu} "  # Define how many core to be used. All cores are now using
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
                f"--limitBAMsortRAM {available_memory} "  # If not set, error can be raised for some cases.
                "--outSAMattributes All XS "
                "--quantMode GeneCounts "
                "--twopassMode Basic "
                "> report_genome_alignment.log"
            ), shell=True)  # The specified command will be executed through the shell.
            raw_bam = os.path.join(job_dir, prefix + "Aligned.sortedByCoord.out.bam")
            
            if not jobbing["is_umi_extracted"]:
                raw_sam = os.path.join(job_dir, prefix + "Aligned.sortedByCoord.out.sam")
                run_and_check(f"cd {job_dir}; {shutil.which('/usr/bin/samtools')} view -h -o {raw_sam} {raw_bam}",
                              shell=True)  # The specified command will be executed through the shell.
                self.jobs[job_id]["process_genome_alignment"] = raw_sam
            else:
                after_sort = os.path.splitext(raw_bam)[0] + '_sorted.bam'
                run_and_check((  # Sort bam file and Index
                    f"cd {job_dir}; "
                    f"{shutil.which('/usr/bin/samtools')} sort {raw_bam} -o {after_sort}; "
                    f"{shutil.which('/usr/bin/samtools')} index {after_sort}"
                ), shell=True)  # The specified command will be executed through the shell.

                # Deduplication of UMI
                output_prefix = "umi-deduplicated"
                print(f"UMI deduplication for {job_id} is now running.")
                run_and_check((  # Run deduplication
                    f"cd {job_dir}; "
                    f"{shutil.which('umi_tools')} dedup "  # Define which umi_tools installation to use
                    f"-I {after_sort} "  # Input
                    f"--output-stats={output_prefix} "
                    f"-S {output_prefix}.bam "
                    f"> {output_prefix}.log "
                    f"> report_{output_prefix}.log"
                ), shell=True)

                run_and_check((  # Convert to sam file
                    f"cd {job_dir}; "
                    f"{shutil.which('/usr/bin/samtools')} view -h "  # Define which samtools installation to use
                    f"-o {output_prefix}.sam {output_prefix}.bam"  # Output & Input
                ), shell=True)  # The specified command will be executed through the shell.

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
        run_and_check((
            f"cd {job_dir}; "  
            f"{shutil.which('bowtie2')} "  # Run Bowtie2 module
            "-D 40 -R 6 -N 0 -L 15 -i S,1,0.50 "  # Alignment effort and sensitivity. It is now very high.
            # Add warning!
            "-I20 -X120 "  # Search only those that has 20-120 nt. 
            "--score-min G,20,5.5 "  # Min score lowered
            "--ma 3 --mp 5,1 "  # ma bonus increased
            "--no-unal "  # To suppress the non-aligned reads
            f"-p{self.cpu} "  # Number of core to use
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
        ), shell=True)  # The specified command will be executed through the shell.

        # Take the transcript sequences into random access memory
        fasta_transcriptome = self.org_db.get_db("cdna")
        transcript_seqs = dict()  # Create a dictionary with transcript id as keys
        with open(fasta_transcriptome, "r") as handle:  # Use previously filtered fasta file
            for record in SeqIO.parse(handle, "fasta"):
                transcript_seqs[record.id] = str(record.seq)  # Sequence as the values

        output_fasta = os.path.join(job_dir, "footprints.fasta")

        # Processing the SAM file
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
            run_and_check((
                f"cd {job_dir}; "  # Change the directory to the index directory
                f"{shutil.which('bowtie2')} "  # Run Bowtie2 module
                f"-p{self.cpu} "  # Number of core to use
                "--no-mixed "  # Do not search for individual pairs if one in a pair does not align.
                # Add warning!
                # Delete below, no need to be biased!
                "-I20 -X120 "  # Default -I=0, -X=500. Since I will disregard X>120 and I<20 in the link-pairing module
                "--time "  # Print the wall-clock time required to load the index files and align the reads.
                "--score-min G,20,6 --ma 4 "  # Allow looser alignment. --ma 3 
                "--local --sensitive-local "  # Apply soft clipping when aligning. Default sensitivity.
                "-k1 "  # We are not interested in best alignment as long as it aligns somewhere in the indexed fasta.
                f"-x {self.index_directories['index_rrna']} "  # Index directory with the base name
                "-q "  # Indicates the inputs are fastq
                f"-1 {temp_paths[0]} "  # Read 1
                f"-2 {temp_paths[1]} "  # Read 2
                f"--un-conc {output_path} "  # Output fastq file, Contains all reads which did not aligned RNAs. 
                # f"--al-conc {os.path.join(job_dir, 'Read%_only_rRNA.fastq')} "  # For testing purposes
                "-S /dev/null "  # Discard alignment sam file /dev/null
                f"2> report_cleanup.log"
            ), shell=True)  # The specified command will be executed through the shell.

            self.jobs[job_id]["process_cleanup"] = [os.path.join(job_dir, j)
                                                    for j in ["Read1_norRNA.fastq", "Read2_norRNA.fastq"]]

        elif jobbing["sequencing_method"] in ["single"]:
            output_path = os.path.join(job_dir, "Read1_norRNA.fastq")
            print(f"rRNA removal for {job_id} is now running.")
            run_and_check((
                f"cd {job_dir}; "  # Change the directory to the index directory
                f"{shutil.which('bowtie2')} "  # Run Bowtie2 module
                f"-p{self.cpu} "  # Number of core to use
                "--time "  # Print the wall-clock time required to load the index files and align the reads.
                "--score-min G,20,6  --ma 4 "  # Allow looser alignment.  --ma 3 
                "--local --sensitive-local "  # Apply soft clipping when aligning. Default sensitivity.
                "-k1 "  # We are not interested in best alignment as long as it aligns somewhere in the indexed fasta.
                f"-x {self.index_directories['index_rrna']} "  # Index directory with the base name
                "-q "  # Indicates the inputs are fastq
                f"{temp_paths} "  # Read 1
                f"--un {output_path} "  # Output fastq file, Contains all reads which did not aligned RNAs. 
                # f"--al {os.path.join(job_dir, 'Read1_only_rRNA.fastq')} "  # For testing purposes
                "-S /dev/null "  # Discard alignment sam file /dev/null
                f"2> report_rnaremove.txt"
            ), shell=True)  # The specified command will be executed through the shell.

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
        assert len(self.jobs) <= self.cpu
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
            pattern1, pattern2 = jobbing["pattern_umi1"], jobbing["pattern_umi2"]
            if not pattern1 or not pattern2:
                return temp_paths
            final_paths = [f"{i}_no-adapt_umi-aware.fastq.gz" for i in ["read1", "read2"]]
            final_paths = [os.path.join(job_dir, i) for i in final_paths]
            print(f"Umitool for {job_id} is now running.")
            run_and_check((
                f"cd {job_dir}; "  # Change the directory to the index directory
                f"{shutil.which('umi_tools')} extract "  # Define which extract installation to use
                "--extract-method=regex "
                f"--bc-pattern='{pattern1}' "  # Barcode pattern
                f"--bc-pattern2='{pattern2}' "  # Barcode pattern for paired reads"
                f"-I {temp_paths[0]} -S {final_paths[0]} "  # Input and output for read 1
                f"--read2-in={temp_paths[1]} --read2-out={final_paths[1]} "  # Input and output for read 2
                f"--log=umi_tools.log "  # Log the results
                f"> 'report_umitool.log'"
            ), shell=True)  # The specified command will be executed through the shell.

            return final_paths

        elif jobbing["sequencing_method"] in ["single"]:
            umitools_pattern = jobbing["pattern_umi"]
            if not umitools_pattern:
                return temp_paths
            final_path = os.path.join(job_dir, "read1_no-adapt_umi-aware.fastq.gz")
            print(f"Umitool for {job_id} is now running.")
            run_and_check((
                f"cd {job_dir}; "  # Change the directory to the index directory
                f"{shutil.which('umi_tools')} extract "  # Define which extract installation to use
                "--extract-method=regex "
                f"--bc-pattern='{umitools_pattern}' "  # Barcode pattern
                f"-I {temp_paths} -S {final_path} "  # Input and output for read 1
                f"--log=umi_tools.log "  # Log the results
                f"> 'report_umitool.log'"
            ), shell=True)  # The specified command will be executed through the shell.

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
            run_and_check((
                    f"cd {job_dir}; "  # Change the directory to the index directory
                    f"{shutil.which('cutadapt')} "  # Define which cutadapt installation to use
                    f"--cores={self.cpu} "  # Define how many core to be used. All cores are now using
                    "--match-read-wildcards "
                    f"-q20 -m23 "  # Settings for quality and minimum length to cut the adapter sequence
                    # todo: explain why: check the photo taken in 12 february
                    # "--discard-untrimmed "
                    # "-O6 "  # Require minimum length overlap between read and adapter for an adapter to be found.
                    f"{read1_adapter} {read2_adapter}".strip() + " "  # Adapter flags. No flanking white space allowed
                    f"-o {temp_paths[0]} -p {temp_paths[1]} "  # Path to output trimmed sequences
                    f"{input_fastq[0]} {input_fastq[1]} "  # Input file
                    f"1> 'report_cutadapt_temp.log'"
            ), shell=True)  # The specified command will be executed through the shell.

            self.jobs[job_id]["process_cutadapt"] = temp_paths

        elif jobbing["sequencing_method"] in ["single"]:
            temp_path = os.path.join(job_dir, "read1_cutadapt_temp.fastq.gz")  # Outputs for the first run
            print(f"Cutadapt for {job_id} is now running.")
            run_and_check((
                f"cd {job_dir}; "  # Change the directory to the index directory
                f"{shutil.which('cutadapt')} "  # Define which cutadapt installation to use
                f"--cores={self.cpu} "  # Define how many core to be used. All cores are now using
                "--match-read-wildcards "
                f"-q20 -m23 "  # Settings for quality and minimum length to cut the adapter sequence
                # If we don't below two, UMI tool definitely malfunction
                "--discard-untrimmed "
                "-O6 "  # Require minimum length overlap between read and adapter for an adapter to be found.
                f"-a {jobbing['adapter']} "  # The adapter flags. No flanking white space allowed
                f"-o {temp_path} "  # Path to output trimmed sequences
                f"{input_fastq} "  # Input file
                f"1> 'report_cutadapt_temp.log'"
            ), shell=True)  # The specified command will be executed through the shell.

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
            run_and_check((
                f"cd {dir_cdna}; "  # Change the directory to the index directory
                f"{shutil.which('bowtie2-build')} "  # Name of the function
                f"--threads {self.cpu} "  # Number of threads to be used.
                f"{temp_cdna_fasta} "  # Input file. -f is to indicate the file is in fasta format
                f"{index_name_cdna} "  # The basename of the index files to write
                f"> report_{self.organism}_cdna.log"
            ), shell=True)  # The specified command will be executed through the shell.
            metadata_index = get_files_metadata(dir_cdna)  # Write the file info
            joblib.dump(metadata_index, os.path.join(dir_cdna, ".metadata.joblib"))

        if not self.create_index_search_metadata(dir_rrna):
            print("Indexing for rRNA is being calculated.")
            temp_rrna_fasta = self.org_db.get_db("rrna")
            run_and_check((
                f"cd {dir_rrna}; "  # Change the directory to the index directory
                f"{shutil.which('bowtie2-build')} "  # Name of the function
                f"--threads {self.cpu} "  # Number of threads to be used.
                f"{temp_rrna_fasta} "  # Input file. -f is to indicate the file is in fasta format
                f"{index_name_rrna} "  # The basename of the index files to write
                f"> report_{self.organism}_rrna.log"
            ), shell=True)  # The specified command will be executed through the shell.
            metadata_index = get_files_metadata(dir_rrna)  # Write the file info
            joblib.dump(metadata_index, os.path.join(dir_rrna, ".metadata.joblib"))

        if not self.create_index_search_metadata(dir_dna):
            # Manual says: "In most cases, the default value of 100 will work as well as the ideal value."
            read_length_minus_1 = 100
            print("Indexing for DNA is being calculated.")
            temp_dna_fasta = self.org_db.get_db("dna")
            temp_gtf_fasta = self.org_db.get_db("gtf")

            # Calculate genomeSAindexNbases
            genome_size = sum([len(rec.seq) for rec in SeqIO.parse(temp_dna_fasta, "fasta")])
            sa_index_parameter = min(14, math.floor(math.log(genome_size, 2) / 2 - 1))
            # For small genomes, this paramter has to be scaled down.
            run_and_check((
                f"cd {dir_dna}; "  # Change the directory to the index directory
                f"{shutil.which('STAR')} "  # Define which star installation to use
                f"--runThreadN {self.cpu} "  # Define how many core to be used. All cores are now using
                "--runMode genomeGenerate "
                f"--genomeDir {dir_dna} "  # Directory to save the files
                f"--genomeSAindexNbases {sa_index_parameter} "
                f"--genomeFastaFiles {temp_dna_fasta} "  # Specifies FASTA file with the genome reference sequences
                f"--sjdbGTFfile {temp_gtf_fasta} "  # Specifies the file with annotated transcripts in the GTF format
                f"--sjdbOverhang {read_length_minus_1} "  # The length of the sequence around the annotated junction
                f"> report_{self.organism}_dna.log"
            ), shell=True)  # The specified command will be executed through the shell.
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

    ID_MAP = {"homo_sapiens": 9606,
              "mus_musculus": 10090,
              "saccharomyces_cerevisiae": 4932,
              "escherichia_coli": 562}

    def __init__(self, organism, ensembl_release, temp_repo_dir):
        assert organism in ["homo_sapiens", "mus_musculus", "saccharomyces_cerevisiae", "escherichia_coli"]
        self.ensembl_release = ensembl_release
        self.organism = organism
        self.organism_id = OrganismDatabase.ID_MAP[self.organism]
        self.temp_repo_dir = temp_repo_dir
        self.repository = create_dir(self.temp_repo_dir, self.organism)

        if self.organism != "escherichia_coli":
            assert 90 <= self.ensembl_release <= 104, "Ensembl Release must be between 90 and 104."
            base_temp = f"ftp://ftp.ensembl.org/pub/release-{ensembl_release}"
        else:
            assert self.ensembl_release == 48, "Ensembl Release must be '48' for escherichia_coli."
            base_temp = f"ftp://ftp.ensemblgenomes.org/pub/release-{ensembl_release}/bacteria"

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
            # Protein Fasta
            pep_temp = "fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
            self.pep = os.path.join(base_temp, pep_temp)
            # Gerp Conservation score
            gerp_temp = "compara/conservation_scores/111_mammals.gerp_conservation_score/" \
                        "gerp_conservation_scores.homo_sapiens.GRCh38.bw"
            self.gerp = os.path.join(base_temp, gerp_temp)
        
        elif self.organism == "mus_musculus":
            # Genome GTF
            gtf_temp = f"gtf/mus_musculus/Mus_musculus.GRCm38.{self.ensembl_release}.chr_patch_hapl_scaff.gtf.gz"
            self.gtf = os.path.join(base_temp, gtf_temp)
            # Genome GFF3
            gff3_temp = f"gff3/mus_musculus/Mus_musculus.GRCm38.{self.ensembl_release}.chr_patch_hapl_scaff.gff3.gz"
            self.gff3 = os.path.join(base_temp, gff3_temp)
            # Genome DNA fasta
            dna_temp = "fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
            self.dna = os.path.join(base_temp, dna_temp)
            # Transcriptome DNA fasta
            cdna_temp = "fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"
            self.cdna = os.path.join(base_temp, cdna_temp)
            # Protein Fasta
            pep_temp = "fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz"
            self.pep = os.path.join(base_temp, pep_temp)
            # Gerp Conservation score
            gerp_temp = "compara/conservation_scores/111_mammals.gerp_conservation_score/" \
                        "gerp_conservation_scores.mus_musculus.GRCm38.bw"
            self.gerp = os.path.join(base_temp, gerp_temp)

        elif self.organism == "saccharomyces_cerevisiae":
            # Genome GTF
            gtf_temp = f"gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.{self.ensembl_release}.gtf.gz"
            self.gtf = os.path.join(base_temp, gtf_temp)
            # Genome GFF3
            gff3_temp = f"gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.{self.ensembl_release}.gff3.gz"
            self.gff3 = os.path.join(base_temp, gff3_temp)
            # Genome DNA fasta
            dna_temp = "fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz" 
            # Note: If the primary assembly file is not present, that indicates there are no haplotype or patch
            # regions, and the 'toplevel' file is equivalent
            self.dna = os.path.join(base_temp, dna_temp)
            # Transcriptome DNA fasta
            cdna_temp = "fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
            self.cdna = os.path.join(base_temp, cdna_temp)
            # Protein Fasta
            pep_temp = "fasta/saccharomyces_cerevisiae/pep/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz"
            self.pep = os.path.join(base_temp, pep_temp)

        elif self.organism == "escherichia_coli":
            # Genome GTF
            gtf_temp = "gtf/bacteria_15_collection/escherichia_coli_bl21_de3_gca_000022665/" \
                       "Escherichia_coli_bl21_de3_gca_000022665.ASM2266v1.48.gtf.gz"
            self.gtf = os.path.join(base_temp, gtf_temp)
            # Genome GFF3
            gff3_temp = "gff3/bacteria_15_collection/escherichia_coli_bl21_de3_gca_000022665/" \
                        "Escherichia_coli_bl21_de3_gca_000022665.ASM2266v1.46.gff3.gz"
            self.gff3 = os.path.join(base_temp, gff3_temp)
            # Genome DNA fasta
            dna_temp = "fasta/bacteria_15_collection/escherichia_coli_bl21_de3_gca_000022665/dna/" \
                       "Escherichia_coli_bl21_de3_gca_000022665.ASM2266v1.dna.nonchromosomal.fa.gz"
            # Note: Nonchromosomal contains DNA that has not been assigned a chromosome
            self.dna = os.path.join(base_temp, dna_temp)
            # Transcriptome DNA fasta
            cdna_temp = "fasta/bacteria_15_collection/escherichia_coli_bl21_de3_gca_000022665/cdna/" \
                        "Escherichia_coli_bl21_de3_gca_000022665.ASM2266v1.cdna.all.fa.gz"
            self.cdna = os.path.join(base_temp, cdna_temp)
            # Protein Fasta
            pep_temp = "fasta/bacteria_15_collection/escherichia_coli_bl21_de3_gca_000022665/pep/" \
                       "Escherichia_coli_bl21_de3_gca_000022665.ASM2266v1.pep.all.fa.gz"
            self.pep = os.path.join(base_temp, pep_temp)

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
            if not os.access(output_path_compressed, os.R_OK) or not os.path.isfile(output_path_compressed):
                print(f"Downloading from the server for {db}:{os.linesep}{db_url}")
                run_and_check(f"cd {self.repository}; curl -L -O --silent {db_url}", shell=True)
            subprocess.run(f"cd {self.repository}; gzip -d -q {output_path_compressed}", shell=True)
        return output_path_uncompressed

    def get_uncompressed_db(self, db):
        db_url = eval(f"self.{db}")
        output_path = os.path.join(self.repository, os.path.basename(db_url))
        if not os.access(output_path, os.R_OK) or not os.path.isfile(output_path):
            print(f"Downloading from the server for {db}:{os.linesep}{db_url}")
            run_and_check(f"cd {self.repository}; curl -L -O --silent {db_url}", shell=True)
        return output_path

    def _download_rna_central(self, db_url):
        output_path = os.path.join(self.repository, os.path.basename(db_url))
        if not os.access(output_path, os.R_OK) or not os.path.isfile(output_path):
            print(f"Downloading from the server for rrna:{os.linesep}{db_url}")
            run_and_check(f"cd {self.repository}; curl -L -O --silent {db_url}", shell=True)
            # The specified command will be executed through the shell.
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


def run_and_check(the_string, *args, **kwargs):
    s = subprocess.run(the_string, *args, **kwargs)
    if s.returncode != 0:
        print(f"{Col.F}Error in the following subprocess:{os.linesep}{the_string}{os.linesep}{Col.E}")


if __name__ == '__main__':
    main()


# End
