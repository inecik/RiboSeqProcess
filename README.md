# RiboSeqAnalysis

Deep sequencing data processing pipeline. 

### Usage

```
>python3 pipeline.py -h

usage: pipeline.py [-h] [-r IDENTIFIER]
                   [-a {homo_sapiens,mus_musculus,saccharomyces_cerevisiae,escherichia_coli}]
                   [-e ENSEMBL_RELEASE] [-c CPU] [-s {3,5}] -f FILEPATH -t
                   TEMP -o OUTPUT

Deep sequencing data processing pipeline.

optional arguments:
  -h, --help            show this help message and exit
  -r IDENTIFIER         the identifier for this run, which will be the
                        directory name for all outputs under provided main
                        output directory. If skipped, the default is a string
                        containing date and time of the program start.
  -a {homo_sapiens,mus_musculus,saccharomyces_cerevisiae,escherichia_coli}
                        organism of interest for the analysis. If skipped, the
                        default value is homo_sapiens.
  -e ENSEMBL_RELEASE    ensembl version to be used. Ignored if the 'organism'
                        does not necessitates Ensembl release information.
                        Default value is 102. For E. coli, it must be 48.
  -c CPU                number of cpu cores to be used. Default value is
                        maximum minus eight.
  -s {3,5}              select 3' or 5' assignment for Julia script.
  -f FILEPATH           path of the txt file which contains the task list.
  -t TEMP               absolute path of the directory to be used for
                        temporary files such as genome indexes.
  -o OUTPUT             absolute path of the directory to be used for output
                        files.
```

Example terminal command:

```
python3 /home/kai/HDD_Kemal/from_raf_computer/Kemal/RiboSeqProcess/pipeline.py -r JaroData -a escherichia_coli -e 48 -f /home/kai/HDD_Kemal/from_raf_computer/Kemal/RiboSeqProcess/example/jaro.txt -t /home/kai/HDD_Kemal/SequencingProcess/Temp/ -o /home/kai/HDD_Kemal/SequencingProcess/Output/
```

Example task list file can be found under `example` folder of this repository.