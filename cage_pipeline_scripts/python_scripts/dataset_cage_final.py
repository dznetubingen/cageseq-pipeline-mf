#!/usr/bin/env python2.7
import os
import sys
import subprocess
import commands
import argparse
import shutil

# print date, to monitor overall timing
subprocess.call("echo 'Start time' && date", shell=True)
cores="1"
## pathto_bin="/home/vsac/pipelines2/cage"
pathto_bin="/data-analysis/home/francescattom/bin/cage_pipeline_scripts"

#####################################################################
# parse command line arguments
#####################################################################
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory", help="specify directory to run the analysis (e.g. /home/user/)", required=True)
parser.add_argument("-r", "--reference", help="specify the genome reference (e.g. hg19)", required=True)
args = parser.parse_args()

#####################################################################
# define some variables used throughout
#####################################################################
# from parsed arguments get directory and reference
DIR = args.directory
reference = args.reference
# once I have the DIR I can defined who the barcodeFile is
barcodeFile = DIR + "/barcodes.txt"

#####################################################################
# define BS genome, STAR genome dir and feature counts dir depending
# on the reference
#####################################################################
if reference == 'hg38':
    BSgenome = 'BSgenome.Hsapiens.UCSC.hg38'
elif reference == 'hg19':
    BSgenome = 'BSgenome.Hsapiens.UCSC.hg19'
elif reference == 'mm10':
    BSgenome = 'BSgenome.Mmusculus.UCSC.mm10'
elif reference == 'danRer10':
    BSgenome = 'BSgenome.Drerio.UCSC.danRer10'
elif reference == 'dm6':
    BSgenome = 'BSgenome.Dmelanogaster.UCSC.dm6'
elif reference == 'ce10':
    BSgenome = 'BSgenome.Celegans.UCSC.ce10'
else:
    print "Sorry, reference " + reference + " not recognised"

# create directory tree where output will be orderly stored
os.mkdir(DIR+"/input")
os.mkdir(DIR+"/logs")
os.mkdir(DIR+"/QC")
os.mkdir(DIR+"/QC/FastQC")
os.mkdir(DIR+"/needless_files")
os.mkdir(DIR+"/plots")
os.mkdir(DIR+"/stats")
os.mkdir(DIR+"/bam_files")
os.mkdir(DIR+"/ctss_files")

# define the initial fastq file, delete it if already existing
fileName = "00_cage.fastq.gz"
if os.path.exists(fileName):
    os.remove(fileName)

#####################################################################
# start the actual processing
#####################################################################
# create single fastq file with standard name
subprocess.call("cat " + DIR + "/*fastq.gz >> " + DIR + "/00_cage.fastq.gz", shell=True)
# run fastqc on the raw data to have an initial quality assessment
subprocess.call("fastqc -t "+cores+" --extract " + DIR + "/00_cage.fastq.gz  > " + DIR + "/fastqc_all.log 2>&1", shell=True)
# extract sequence length and initial number of reads from fastqc report
seq_length = commands.getoutput("cat " + DIR + "00_cage_fastqc/fastqc_data.txt | grep Sequence | grep length | cut -f 2")
totalReads = commands.getoutput("cat " + DIR + "00_cage_fastqc/fastqc_data.txt | grep 'Total Sequences' | cut -f 2")
# write the total number of reads in the file that will be used to make the plots
value = ["totalReads", str(totalReads), "\n"]
s="\t".join(value)
f = open(DIR+"/stats00_barcode_splitting.txt", 'w')
f.write(s)
f.close()

# store the initial fastqc report as reference
shutil.move(DIR+"/00_cage_fastqc.html", DIR+"/QC/fastqc_report_raw_data.html")
os.remove(DIR+"/00_cage_fastqc.zip")
shutil.rmtree(DIR+"/00_cage_fastqc")

#####################################################################
# iterate through barcode/sample pairs to submit one job for each
# barcode
#####################################################################
with open(barcodeFile) as input_file:
    for i, line in enumerate(input_file):
        line = line.rstrip('\n\r')
        (sampleId, barcode) = line.split("\t")
        subprocess.call("sbatch --nodelist=tu-svr-003 --no-requeue -e " +DIR+"/"+sampleId+ ".std.err -o "+
                        DIR+"/"+sampleId + ".std.out "+pathto_bin+"/python_scripts/sample_cage_final.py -b " +barcode+ " -s " + sampleId + 
                        " -d " +DIR+ " -r " + reference + " -l " + seq_length + " -bs "+BSgenome, shell=True)
