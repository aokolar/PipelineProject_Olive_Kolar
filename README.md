# PipelineProject_Olive_Kolar

## Running this pipeline

### Dependencies:
Python modules:
- os
- shutil
- csv
- statistics
- subprocess
- Biopython SeqIO

Tools:
- kallisto
- sleuth (R package)
- Bowtie2
- SPAdes
- BLAST+

### Other necessary information:

You must download the R script "pipeline_sleuth.R" included in this repo in order to run sleuth. Upload this file into your workspace and keep it in your home directory. Do not put it in a folder.


### Input Data:

#### If you want full-size genomes (this pipeline will take quite some time to run):

Execute the following commands in your terminal, one at a time.

    $ wget https://www.ncbi.nlm.nih.gov/sra/SRX2896360
    
    $ wget https://www.ncbi.nlm.nih.gov/sra/SRX2896363
    
    $ wget https://www.ncbi.nlm.nih.gov/sra/SRX2896374
    
    $ wget https://www.ncbi.nlm.nih.gov/sra/SRX2896375
    
    $ fasterq-dump SRX2896360
    
    $ fasterq-dump SRX2896363
    
    $ fasterq-dump SRX2896374
    
    $ fasterq-dump SRX2896375

This will download and unpack all necessary files. From here, just execute the Python pipeline script and it will do everything else for you. If you want to execute the script through the command line, use this command:

    python pipeline_wrapper.py



#### For smaller sample-size genomes (this pipeline will run in under 1 minute):

Download all of the files named "SRX_______.fastq" (there are 8 of them) in this repo and upload them into your workspace. Keep them in your home directory, and do not put them in a folder.

From there, just execute the Python pipeline script and it will do everything else for you. If you want to execute the script through the command line, use this command:

    python pipeline_wrapper.py



### Output Information:

Details about what happened during the pipeline's run can be found in your "project_files" directory in the file "PipelineProject.log". A completed log from a run of this tool with full-size genomes is included in this repository.
