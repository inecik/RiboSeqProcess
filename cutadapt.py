"""
Cutadapt for paired end sequencing
"""


from subprocess import run
from multiprocessing import cpu_count
from shutil import which

def cutadapt(inputs, outputs, trim_forward, trim_reverse, )
    trim_flanking = [-2, 5]




    cut_adapt_command = (
        f"cutadapt --cores={n_core} -u 2 -U 5 -o "  # settings for cutadapt
        f"'{trim_outputs[0]}' -p '{trim_outputs[1]}' '{excel_path_s[0]}' '{excel_path_s[1]}' "  # outputs and inputs
        f"1> '{trim_output_dir}/{excel_rep}_Cutadapt_report.txt"  # Export the logs
    )
subprocess.run(cut_adapt_command, shell=True)


if __name__ == "__main__":
    my_function()
