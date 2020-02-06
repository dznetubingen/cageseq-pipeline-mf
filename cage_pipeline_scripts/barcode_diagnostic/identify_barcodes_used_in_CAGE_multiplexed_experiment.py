#!/usr/bin/env python2.7
import os
import sys
import subprocess
import commands
import argparse
import shutil
import glob

# print date, to monitor overall timing
subprocess.call("echo 'Start time' && date", shell=True)

## set the number of cores to be used for each skewer run (note that
## they are run sequentially, so you can use a lot of CPUs if nothing
## else is running)
cores="10"
## this is where all possible barcodes are stored
barcodeFile="/data-analysis/home/francescattom/bin/cage_pipeline_scripts/barcode_diagnostic/all_possible_barcodes.txt"

#####################################################################
# parse command line arguments
#####################################################################
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory", help="specify directory to run the analysis (e.g. /home/user/)", required=True)
args = parser.parse_args()

# from parsed arguments get directory and reference
DIR = args.directory

## # create directory tree where output will be orderly stored
os.mkdir(DIR+"/input")
os.mkdir(DIR+"/QC")
os.mkdir(DIR+"/QC/FastQC")
os.mkdir(DIR+"/needless_files")

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
# extract sequence length and store it to run the barcode_splitting
seq_length = commands.getoutput("cat " + DIR + "00_cage_fastqc/fastqc_data.txt | grep Sequence | grep length | cut -f 2")
maxL1 = str(int(seq_length) - 3)
minL1 = maxL1
# extract initial number of reads from fastqc report to plot the barcode splitting plot
totalReads = commands.getoutput("cat " + DIR + "00_cage_fastqc/fastqc_data.txt | grep 'Total Sequences' | cut -f 2")
# write the total number of reads in the file that will be used to make the plots
value = ["totalReads", str(totalReads), "\n"]
s="\t".join(value)
f = open(DIR+"/stats00_barcode_splitting.txt", 'w')
f.write(s)
f.close()

# store the initial fastqc report as reference
shutil.move(DIR+"/00_cage_fastqc.html", DIR+"/QC/fastqc_report_raw_data.html")
# remove the rest
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
        subprocess.call("module load skewer; skewer -z -t "+cores+" -m head -x " + barcode + " -r 0 -d 0 -k 3 -l " + 
                        minL1 + " -L " + maxL1 + " " + DIR + "/00_cage.fastq.gz --quiet -o "+
                        DIR+"/01_" + sampleId + " > "+DIR+"/log_all_skewer_std_out_" + sampleId + 
                        " 2>&1", shell=True)
        subprocess.call("fastqc --extract -t "+cores+" "+DIR+"/01_" + sampleId + 
                        "-trimmed.fastq.gz >> "+DIR+"/fastqc_all.log 2>&1", shell=True)
        readsInitial = commands.getoutput("cat " + DIR + "01_" + sampleId + 
                                          "-trimmed_fastqc/fastqc_data.txt | grep 'Total Sequences' | cut -f 2")
        value = [barcode, str(readsInitial), "\n"]
        s="\t".join(value)
        f = open(DIR+"/stats00_barcode_splitting.txt", 'a')
        f.write(s)
        f.close()

        # clean up the fastqc run
        shutil.move(DIR+"/01_" + sampleId + "-trimmed_fastqc.html", DIR+"/QC/FastQC/")
        shutil.move(DIR+"/01_" + sampleId + "-trimmed.fastq.gz", DIR+"/needless_files/")
        shutil.move(DIR+"/01_" + sampleId + "-trimmed.log", DIR+"/needless_files/")
        shutil.move(DIR+"/log_all_skewer_std_out_" + sampleId, DIR+"/needless_files/")
        os.remove(DIR+"/01_" + sampleId + "-trimmed_fastqc.zip")
        shutil.rmtree(DIR+"/01_" + sampleId + "-trimmed_fastqc")
        # inform user which barcode we are at
        print(sampleId)
        sys.stdout.flush()

## run the R script that makes the barcode splitting plot
subprocess.call("/data-analysis/home/francescattom/bin/cage_pipeline_scripts/R_scripts/r00_barcode_splitting.R "+
                DIR+"/stats00_barcode_splitting.txt "+DIR+"/plot00_barcode_splitting_barplot.png", shell=True)

## final cleanup
os.remove(DIR+"/00_cage.fastq.gz")
os.remove(DIR+"/fastqc_all.log")

## these files need to be cancelled in the "real" version
for filename in glob.glob(DIR+"/*fastq.gz"):
    shutil.move(filename, DIR+"/input/")

# print date, to monitor overall timing
subprocess.call("echo 'End time' && date", shell=True)


