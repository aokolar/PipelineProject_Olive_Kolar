import os
import shutil
import csv
import statistics
import subprocess
from Bio import SeqIO


#create folder to hold files for this pipeline
make_folder = "mkdir project_files"


######## KALLISTO ########


#command to download reference transcriptome CDS to build kallisto index from
download_genome = "datasets download virus genome accession NC_006273.2 --include cds --filename virus_genome.zip"
#unzip files from this download and put them in project_files
unzip = "unzip virus_genome.zip -d project_files"

#run these commands in the command line
os.system(make_folder)
os.system(download_genome)
os.system(unzip)


#name for the index that will be created by kallisto
index_for_kallisto = "index.idx"
#reference transcriptome to build index from
ref_transcriptome = "cds.fna"

#navigate to the folder "data" which holds the CDS file
os.chdir("project_files/ncbi_dataset/data")

#creating the kallisto index-building command
kallisto_index = "kallisto index -i " +index_for_kallisto+ " " +ref_transcriptome

#run kallisto index command
os.system(kallisto_index)


#open the FASTA file containing the CDS from the reference transcriptome
fastafile = open("cds.fna")
n = 0
#count number of CDS in file
for line in fastafile:
    if line.startswith(">"):
        n += 1
fastafile.close()


#navigate back to project_files
os.chdir("..")
os.chdir("..")


#create log file
f = open("PipelineProject.log", "w")
f.close()

#write number of CDS in reference transcriptome to log file
f = open("PipelineProject.log", "a")
f.write(f"The HCMV genome (NC_006273.2) has {n} CDS.")
f.close()


#navigate back to home directory
os.chdir("..") 



#move index file from nested "data" folder back to just the "project_files" folder so it's easier to access
shutil.move("project_files/ncbi_dataset/data/index.idx", "project_files/index.idx")
#move all paired-end reads to "project_files" folder so everything is in the same place
shutil.move("SRX2896360_1.fastq", "project_files/SRX2896360_1.fastq")
shutil.move("SRX2896360_2.fastq", "project_files/SRX2896360_2.fastq")
shutil.move("SRX2896363_1.fastq", "project_files/SRX2896363_1.fastq")
shutil.move("SRX2896363_2.fastq", "project_files/SRX2896363_2.fastq")
shutil.move("SRX2896374_1.fastq", "project_files/SRX2896374_1.fastq")
shutil.move("SRX2896374_2.fastq", "project_files/SRX2896374_2.fastq")
shutil.move("SRX2896375_1.fastq", "project_files/SRX2896375_1.fastq")
shutil.move("SRX2896375_2.fastq", "project_files/SRX2896375_2.fastq")


#navigate into project_files folder
os.chdir("project_files")

#create a folder within project_files called "kallisto_results" to hold kallisto outputs
make_kallisto_res_folder = "mkdir kallisto_results"
os.system(make_kallisto_res_folder)


#create kallisto commands for each of the SRX paired-end reads and output results to kallisto_results -> SRX_____
kallisto_quant_1 = "kallisto quant -i index.idx -o kallisto_results/SRX2896360 -b 30 -t 4 SRX2896360_1.fastq SRX2896360_2.fastq"
kallisto_quant_2 = "kallisto quant -i index.idx -o kallisto_results/SRX2896363 -b 30 -t 4 SRX2896363_1.fastq SRX2896363_2.fastq"
kallisto_quant_3 = "kallisto quant -i index.idx -o kallisto_results/SRX2896374 -b 30 -t 4 SRX2896374_1.fastq SRX2896374_2.fastq"
kallisto_quant_4 = "kallisto quant -i index.idx -o kallisto_results/SRX2896375 -b 30 -t 4 SRX2896375_1.fastq SRX2896375_2.fastq"

#run kallisto commands
os.system(kallisto_quant_1)
os.system(kallisto_quant_2)
os.system(kallisto_quant_3)
os.system(kallisto_quant_4)



#create empty txt file that will be used as sleuth input later
f = open("sleuth_table.txt", "w")
f.close()


#now we look at the kallisto results for each sample and collect the TPM values

#navigate to first folder of kallisto results
os.chdir("kallisto_results/SRX2896360")

tpms_1 = [] #list to hold all TPM values in abundance file
with open("abundance.tsv") as file:
    rd = csv.reader(file, delimiter="\t")
    next(file) #skip first line
    for row in rd:
        #add TPM values to list
        tpms_1.append(float(row[4]))

#find the min, median, mean, and max values from the TPM values in the file
min_1 = min(tpms_1)
med_1 = statistics.median(tpms_1)
mean_1 = "%.5f"%(sum(tpms_1) / len(tpms_1))
max_1 = max(tpms_1)


#navigate back to project_files (move up 2 levels)
os.chdir("..")
os.chdir("..")


#navigate to second folder of kallisto results
os.chdir("kallisto_results/SRX2896363")

tpms_2 = [] #list to hold all TPM values in abundance file
with open("abundance.tsv") as file:
    rd = csv.reader(file, delimiter="\t")
    next(file) #skip first line
    for row in rd:
        #add TPM values to list
        tpms_2.append(float(row[4]))

#find the min, median, mean, and max values from the TPM values in the file
min_2 = min(tpms_2)
med_2 = statistics.median(tpms_2)
mean_2 = "%.5f"%(sum(tpms_2) / len(tpms_2))
max_2 = max(tpms_2)


#navigate back to project_files (move up 2 levels)
os.chdir("..")
os.chdir("..")


#navigate to third folder of kallisto results
os.chdir("kallisto_results/SRX2896374")

tpms_3 = [] #list to hold all TPM values in abundance file
with open("abundance.tsv") as file:
    rd = csv.reader(file, delimiter="\t")
    next(file) #skip first line
    for row in rd:
        #add TPM values to list
        tpms_3.append(float(row[4]))

#find the min, median, mean, and max values from the TPM values in the file
min_3 = min(tpms_3)
med_3 = statistics.median(tpms_3)
mean_3 = "%.5f"%(sum(tpms_3) / len(tpms_3))
max_3 = max(tpms_3)


#navigate back to project_files (move up 2 levels)
os.chdir("..")
os.chdir("..")


#navigate to fourth folder of kallisto results
os.chdir("kallisto_results/SRX2896375")

tpms_4 = [] #list to hold all TPM values in abundance file
with open("abundance.tsv") as file:
    rd = csv.reader(file, delimiter="\t")
    next(file) #skip first line
    for row in rd:
        #add TPM values to list
        tpms_4.append(float(row[4]))

#find the min, median, mean, and max values from the TPM values in the file
min_4 = min(tpms_4)
med_4 = statistics.median(tpms_4)
mean_4 = "%.5f"%(sum(tpms_4) / len(tpms_4))
max_4 = max(tpms_4)


#navigate back to project_files (move up 2 levels)
os.chdir("..")
os.chdir("..")


# write TPM stats for each sample to project log file
f = open("PipelineProject.log", "a")
f.write("\n\nKallisto Output:\n")
f.write("sample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n")
f.write(f"SRX2896360\t2dpi\t{min_1}\t{med_1}\t{mean_1}\t{max_1}\n")
f.write(f"SRX2896363\t6dpi\t{min_2}\t{med_2}\t{mean_2}\t{max_2}\n")
f.write(f"SRX2896374\t2dpi\t{min_3}\t{med_3}\t{mean_3}\t{max_3}\n")
f.write(f"SRX2896375\t6dpi\t{min_4}\t{med_4}\t{mean_4}\t{max_3}\n")
f.close()

#write kallisto output info to sleuth table file
f = open("sleuth_table.txt", "a")
f.write("sample condition path\n")
f.write("SRX2896360 2dp1 kallisto_results/SRX2896360\n")
f.write("SRX2896363 6dp1 kallisto_results/SRX2896363\n")
f.write("SRX2896374 2dp1 kallisto_results/SRX2896374\n")
f.write("SRX2896375 6dp1 kallisto_results/SRX2896375\n")
f.close()


####### SLEUTH #######


#navigate to home directory
os.chdir("..")
#move sleuth R script into project_files
shutil.move("pipeline_sleuth.R", "project_files/pipeline_sleuth.R")
#navigate back into project_files
os.chdir("project_files")


#call sleuth R script, run with input from sleuth_table, output results as "fdr05_results" in project_files folder
subprocess.call("Rscript pipeline_sleuth.R", shell=True)

#output sleuth results to log file
f = open("PipelineProject.log", "a")
f.write("\nSleuth Output:\n")
#pull relevant results from sleuth output file and write to log file
with open("fdr05_results.txt") as file:
    rd = csv.reader(file, delimiter=" ")
    for row in rd:
        f.write(f"{row[0]}\t{row[3]}\t{row[1]}\t{row[3]}\n")
f.close()


#navigate back to home directory
os.chdir("..")

#move reference genome file to project_files directory so it is accessible
shutil.move("project_files/ncbi_dataset/data/cds.fna", "project_files/cds.fna")

#navigate to folder containing HCMV reference genome
os.chdir("project_files")



####### BOWTIE2 #######


#command to build index from reference genome using bowtie2
bowtie_build = "bowtie2-build cds.fna HCMV"
#run command to build index
os.system(bowtie_build)


#commands to run bowtie2 for each set of paired-end reads. will output the reads that mapped to the reference genome
run_bowtie_1 = "bowtie2 -x HCMV -1 SRX2896360_1.fastq -2 SRX2896360_2.fastq -S HCMVmap1.sam --al-conc HCMVmapped1.fastq"
run_bowtie_2 = "bowtie2 -x HCMV -1 SRX2896363_1.fastq -2 SRX2896363_2.fastq -S HCMVmap2.sam --al-conc HCMVmapped2.fastq"
run_bowtie_3 = "bowtie2 -x HCMV -1 SRX2896374_1.fastq -2 SRX2896374_2.fastq -S HCMVmap3.sam --al-conc HCMVmapped3.fastq"
run_bowtie_4 = "bowtie2 -x HCMV -1 SRX2896375_1.fastq -2 SRX2896375_2.fastq -S HCMVmap4.sam --al-conc HCMVmapped4.fastq"
os.system(run_bowtie_1)
os.system(run_bowtie_2)
os.system(run_bowtie_3)
os.system(run_bowtie_4)


#count numbers of total reads in each SAM file
f = open("HCMVmap1.sam")
lines = f.readlines()
lastline = lines[-1] #take last line from file
items = lastline.split() #split apart the items in the line
index = items[0] #take the first item from the line
accession, reads = index.split(".") #split apart the SRX accession and the # of reads
sam1_count = reads #keep the number of reads
f.close()

#repeat for other 3 SAM files
f = open("HCMVmap2.sam")
lines = f.readlines()
lastline = lines[-1]
items = lastline.split()
index = items[0]
accession, reads = index.split(".")
sam2_count = reads
f.close()

f = open("HCMVmap3.sam")
lines = f.readlines()
lastline = lines[-1]
items = lastline.split()
index = items[0]
accession, reads = index.split(".")
sam3_count = reads
f.close()

f = open("HCMVmap4.sam")
lines = f.readlines()
lastline = lines[-1]
items = lastline.split()
index = items[0]
accession, reads = index.split(".")
sam4_count = reads
f.close()


#now count number of reads for each sample after bowtie filtering
f = open("HCMVmapped1.1.fastq")
count=0
lines = f.readlines()
for line in lines:
    if line[0]=="@": #count each unique read in the mapped file
        count+=1
map1_count = count #save # of reads for later
f.close()

#repeat for other filtered read files
f = open("HCMVmapped2.1.fastq")
count=0
lines = f.readlines()
for line in lines:
    if line[0]=="@":
        count+=1
map2_count = count
f.close()

f = open("HCMVmapped3.1.fastq")
count=0
lines = f.readlines()
for line in lines:
    if line[0]=="@":
        count+=1
map3_count = count
f.close()

f = open("HCMVmapped4.1.fastq")
count=0
lines = f.readlines()
for line in lines:
    if line[0]=="@":
        count+=1
map4_count = count
f.close()


#add details of read counts to log file
f = open("PipelineProject.log", "a")
f.write("\nBowtie2 Output:\n")
f.write(f"Donor 1 (2dpi) had {sam1_count} read pairs before Bowtie2 filtering and {map1_count} read pairs after.\n")
f.write(f"Donor 1 (6dpi) had {sam2_count} read pairs before Bowtie2 filtering and {map2_count} read pairs after.\n")
f.write(f"Donor 3 (2dpi) had {sam3_count} read pairs before Bowtie2 filtering and {map3_count} read pairs after.\n")
f.write(f"Donor 3 (6dpi) had {sam4_count} read pairs before Bowtie2 filtering and {map4_count} read pairs after.\n")
f.close()




####### SPADES #######


#spades input commands for each donor. uses 2 sets of paired-end reads for each (one pair for each condition), 2 threads, and a k-mer size of 77
spades_input1 = "spades.py --pe-1 1 HCMVmapped1.1.fastq --pe-2 1 HCMVmapped1.2.fastq --pe-1 2 HCMVmapped2.1.fastq --pe-2 2 HCMVmapped2.2.fastq -t 2 -k 77 -o donor1_assembly"
spades_input2 = "spades.py --pe-1 1 HCMVmapped3.1.fastq --pe-2 1 HCMVmapped3.2.fastq --pe-1 2 HCMVmapped4.1.fastq --pe-2 2 HCMVmapped4.2.fastq -t 2 -k 77 -o donor3_assembly"
#run commands
os.system(spades_input1)
os.system(spades_input2)

#write commands used to log file
f = open("PipelineProject.log", "a")
f.write("\nSPAdes Input Commands Used:\n")
f.write(f"{spades_input1}\n")
f.write(f"{spades_input2}\n")
f.close()



####### BLAST+ #######



# find the longest contig in each of the 2 spades output assemblies
os.chdir("donor1_assembly")
max_len_1 = 0
max_description_1 = ""
for record in SeqIO.parse("contigs.fasta", "fasta"): #parse through contigs
    if len(record.seq) > max_len_1:
        max_len_1 = len(record.seq) #keep the longest contig

#write longest contig to a new file called blast1input.fasta
f = open("blast1input.fasta", "w")
for record in SeqIO.parse("contigs.fasta", "fasta"):
    if len(record.seq) == max_len_1:
        f.write(str(record.seq))
f.close()


#navigate back to project_files
os.chdir("..")


#repeat for donor 3 assembly
os.chdir("donor3_assembly")
max_len_3 = 0
max_description_3 = ""
for record in SeqIO.parse("contigs.fasta", "fasta"):
    if len(record.seq) > max_len_3:
        max_len_3 = len(record.seq)

#write longest contig to a new file called blast3input.fasta
f = open("blast3input.fasta", "w")
for record in SeqIO.parse("contigs.fasta", "fasta"):
    if len(record.seq) == max_len_3:
        f.write(str(record.seq))
f.close()
    

#navigate to home directory
os.chdir("..")
os.chdir("..")

#move blast input files to project_files so they are accessible
shutil.move("project_files/donor1_assembly/blast1input.fasta", "project_files/blast1input.fasta")
shutil.move("project_files/donor3_assembly/blast3input.fasta", "project_files/blast3input.fasta")

#navigate back to project_files
os.chdir("project_files")


# remove downloaded files from the old ncbi dataset, they are not needed anymore and will cause hangups when unzipping a new ncbi dataset
os.system("rm -rf ncbi_dataset")
os.system("rm md5sum.txt")
os.system("rm README.md")


#get dataset of Betaherpesvirinae genomes to BLAST against
get_dataset = "datasets download virus genome taxon Betaherpesvirinae --refseq --include genome --filename taxon_dataset.zip"
os.system(get_dataset)
#unzip the dataset. will be put in the "ncbi_dataset/data" directory
os.system("unzip taxon_dataset.zip")

#navigate to home directory
os.chdir("..")

#move dataset file to project_files for easier access
shutil.move("project_files/ncbi_dataset/data/genomic.fna", "project_files/genomic.fna")

#go back to project_files
os.chdir("project_files")

#command to build a database from these genomes, which will be used to BLAST against
make_database = "makeblastdb -in genomic.fna -out betaherpesvirinae -title betaherpesvirinae -dbtype nucl"
os.system(make_database)

#BLAST commands for both assemblies. collects necessary values to put in log file
blast_1 = "blastn -query blast1input.fasta -db betaherpesvirinae -out blast1_output.tsv -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle'"
blast_3 = "blastn -query blast3input.fasta -db betaherpesvirinae -out blast3_output.tsv -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle'"
os.system(blast_1)
os.system(blast_3)


#write outputs to log file
f = open("PipelineProject.log", "a")
f.write("\nBLAST Results:\n")
f.write("Donor 1:\n")
f.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
with open("blast1_output.tsv") as file:
    lines = file.readlines() #access lines in BLAST output (data is already tab-delimited)
    for line in lines:
        f.write(line)
f.write("Donor 3:\n")
f.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
with open("blast3_output.tsv") as file:
    lines = file.readlines()
    for line in lines:
        f.write(line)
f.close()