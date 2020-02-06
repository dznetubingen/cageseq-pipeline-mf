#!/usr/bin/env python2.7

import os
import sys
import subprocess
import commands
import argparse
import shutil
import glob

# parse command line options, currently passed on from p01 script,
# which takes care of creating a large merged fastq file and create
# barcode-oriented jobs
parser = argparse.ArgumentParser()
parser.add_argument("-b", "--barcode", help="specify the barcode for this library (e.g. ATG)", required=True)
parser.add_argument("-s", "--sampleId", help="specify the sample ID of this library (e.g. 10196_caudate)", required=True)
parser.add_argument("-d", "--directory", help="specify directory to run the analysis (e.g. /home/user)", required=True)
parser.add_argument("-r", "--reference", help="specify the genome reference (e.g. hg19))", required=True)
parser.add_argument("-l", "--sequencelength", help="specify the length of the reads (e.g. 51))", required=True)
parser.add_argument("-bs", "--BSgenome", help="", required=True)
args = parser.parse_args()

# assign information read from command line to variables
barcode = args.barcode
sampleId = args.sampleId
DIR = args.directory
reference = args.reference
BSgenome=args.BSgenome
length = int(args.sequencelength)

number_samples= sum(1 for line in open(DIR + "/barcodes.txt"))

genome_dir = "/data-analysis/resources/cage_pipeline/" + reference + "/star"
featureCountsDir = "/data-analysis/resources/cage_pipeline/" + reference + "/featureCounts_files"
pathto_bin="/data-analysis/home/francescattom/bin/cage_pipeline_scripts"

## genome_dir = "/home/vsac/database2/" + reference + "/star"
## featureCountsDir = "/home/vsac/database2/"+reference

#####################################################################
# define sequence length related parameters depending on sequence
# length (passed as command line argument); the first two constraints
# (maxL1, minL1, maxL2 and minL2) are necessary to oblige skewer to
# perform an "anchored 5' trimming" (see cutadapt documentation for an
# explanation of what this means). The last constraint should be
# discussed a bit more maybe: for the moment I leave the minimum
# length to 15bp. The max length at 45 doesn't really make a
# difference since more trimming was done in previous steps in any
# case
#####################################################################
maxL1 = str(length - 3)
minL1 = maxL1
maxL2 = str(length - 9)
minL2 = maxL2
maxL3 = str(45)
minL3 = str(15)
cores="1"

#####################################################################
# first skewer call to remove the three bp barcodes at the 5' end and
# associated fastqc run
#####################################################################
subprocess.call("module load skewer; skewer -z -t "+cores+" -m head -x " + barcode + " -r 0 -d 0 -k 3 -l " + minL1 + " -L " 
                + maxL1 + " " + DIR + "/00_cage.fastq.gz --quiet -o "+DIR+"/01_" + sampleId + 
                " > "+DIR+"/log_all_skewer_std_out_" + sampleId + " 2>&1", shell=True)
subprocess.call("fastqc --extract -t "+cores+" "+DIR+"/01_" + sampleId + "-trimmed.fastq.gz >> "+DIR+"/fastqc_all.log 2>&1", shell=True)

#####################################################################
# collect information to create barcode splitting stats. This is a
# crucial QC information for the demultiplexing, to make sure all
# barcodes provided are correct
#####################################################################
readsInitial = commands.getoutput("cat " + DIR + "01_" + sampleId + "-trimmed_fastqc/fastqc_data.txt | grep 'Total Sequences' | cut -f 2")
value = [sampleId, str(readsInitial), "\n"]
s="\t".join(value)
f = open(DIR+"/stats00_barcode_splitting.txt", 'a')
f.write(s)
f.close()

# clean up the fastqc run
shutil.move(DIR+"/01_" + sampleId + "-trimmed_fastqc.html", DIR+"/QC/FastQC/")
os.remove(DIR+"/01_" + sampleId + "-trimmed_fastqc.zip")
shutil.rmtree(DIR+"/01_" + sampleId + "-trimmed_fastqc")

#####################################################################
# second skewer call to remove 5' adaptor (CAGCAG sequence) + fastqc
# -k 6 requires at least a 6bp overlap to trim
# -r 0 specifies that no error is tolerated
# -d 0 specifies that no indel error is tolerated
# -l, -L: specify min and max length accepted after trimming, 
# necessary to impose that a read is kep only if exactly the sequence 
# "CAGCAG" is trimmed out
#####################################################################
subprocess.call("module load skewer; skewer -z -t "+cores+" -m head -x CAGCAG -k 6 -d 0 -r 0 -l "+minL2+" -L "+maxL2+" " + 
                DIR + "/01_" + sampleId + "-trimmed.fastq.gz --quiet -o "+DIR+"/02_" + sampleId  + 
                " >> "+DIR+"/log_all_skewer_std_out_" + sampleId + " 2>&1", shell=True)
subprocess.call("fastqc -t "+cores+" "+DIR+"/02_" + sampleId + "-trimmed.fastq.gz >> "+DIR+"/fastqc_all.log 2>&1", shell=True)
os.remove(DIR+"/02_" + sampleId + "-trimmed_fastqc.zip")
shutil.move(DIR+"/02_" + sampleId + "-trimmed_fastqc.html", DIR+"/QC/FastQC/")

#####################################################################
# third skewer call, to remove 3'end adapter + fastqc (from which I
# save the sequence length distribution)
#####################################################################
subprocess.call("module load skewer; skewer -t "+cores+" -x TCGTATGCCGTCTTC -l "+minL3+" -L "+maxL3+" "+DIR+"/02_" + sampleId + 
                "-trimmed.fastq.gz --quiet -o "+DIR+"/03_finalTrim_" + sampleId + " >> "+DIR+"/log_all_skewer_std_out_" + sampleId + " 2>&1", shell=True)
subprocess.call("fastqc --extract -t "+cores+" "+DIR+"/03_finalTrim_" + sampleId + "-trimmed.fastq >> "+DIR+"/fastqc_all.log 2>&1", shell=True)

# save the number of reads left after the trimming for QC plot
readsPostTrimming = commands.getoutput("cat " + DIR + "/03_finalTrim_" + sampleId + "-trimmed_fastqc/fastqc_data.txt | grep 'Total Sequences' | cut -f 2")

shutil.move(DIR+"/03_finalTrim_" + sampleId + "-trimmed_fastqc/Images/sequence_length_distribution.png",
            DIR+"/plots/plot03_sequence_length_distribution_after_trimming_"+sampleId+".png")
shutil.move(DIR+"/03_finalTrim_" + sampleId + "-trimmed_fastqc.html", DIR+"/QC/FastQC/")
shutil.rmtree(DIR+"/03_finalTrim_" + sampleId + "-trimmed_fastqc")
os.remove(DIR+"/03_finalTrim_" + sampleId + "-trimmed_fastqc.zip")


###############################################################
## tagdust step
###############################################################
# Note: the artifact inputs for tagdust need to be tailored based on
# the barcode. In particular in the artifact file there is one "nnn"
# string that needs to be replaced by the barcode and an "mmm" string
# that needs to be replaced by the reverse of the barcode. Since there
# are only 10 possible barcodes that are used in CAGE library
# preparation, I hard coded the subsitutions rather than doing
# something fancier
barcode_substitutions = {"ACA":"TGT", "ACG":"TGC", "ACT":"TGA", "AGA":"TCT", "AGC":"TCG", 
                         "AGT":"TCA", "ATC":"TAG", "ATG":"TAC", "CTT":"GAA", "GAC":"CTG", 
                         "GAT":"CTA", "GCC":"CGG", "GTA":"CAT", "TAG":"ATC", "TGG":"ACC"}
# open artifact template, we assume it has been copied in the
# directory where we have input file/s and barcodes
f_in = open(DIR + "/CAGE_artefacts_201507.txt", 'r')
in_artefacts = f_in.read()
f_in.close()
# make the substitution
out_artefacts = in_artefacts.replace("nnn", barcode)
out_artefacts = out_artefacts.replace("mmm", barcode_substitutions[barcode])
# write the tailored artifact file
f_out = open(DIR + "/barcode_" + barcode + "_CAGE_artefacts_201507.txt",'w')
f_out.write(out_artefacts)
f_out.close()
# call tagdust with the tailored command
subprocess.call(pathto_bin+"/tagdust "+DIR+"/barcode_" + barcode + "_CAGE_artefacts_201507.txt "+
                DIR + "/03_finalTrim_" + sampleId + "-trimmed.fastq -a " + 
                DIR + "/tagdust_removed_" + sampleId + ".fastq -o " + DIR + "/04_tagdust_retained_" 
                + sampleId + ".fastq > " + DIR + "/log_tagdust_"+ sampleId +" 2>&1", shell = True)
subprocess.call("fastqc --extract -t "+cores+" "+DIR+"/04_tagdust_retained_" + sampleId + ".fastq >> "+DIR+"/fastqc_all.log 2>&1", shell=True)
shutil.move(DIR+"/04_tagdust_retained_"+sampleId+"_fastqc.html", DIR+"/QC/FastQC/")
shutil.move(DIR+"/04_tagdust_retained_"+sampleId+"_fastqc/Images/per_base_sequence_content.png", DIR+"/plots/plot04_per_base_sequence_content_"+sampleId+".png")
os.remove(DIR+"/04_tagdust_retained_"+sampleId+"_fastqc.zip")
shutil.rmtree(DIR+"/04_tagdust_retained_"+sampleId+"_fastqc")

################################################################################
# STAR call for the mapping and fastqc (from which we save the full
# report and the duplication levels plot, that gives information
# related to library redundancy)
################################################################################
subprocess.call("STAR --outFileNamePrefix "+DIR+"/star_" + sampleId + "_ --genomeDir "+ genome_dir + " --outFilterMismatchNoverLmax 0.08"
                " --alignIntronMax 1 --runThreadN "+cores+" --outFilterMatchNmin 15 --outFilterScoreMinOverLread 0"
                " --outFilterMatchNminOverLread 0  --readFilesIn "+DIR+"04_tagdust_retained_" + sampleId + ".fastq"
                " --limitBAMsortRAM 10000000000 --genomeLoad LoadAndRemove --alignEndsType EndToEnd"
                " --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outFilterMultimapNmax 1 > "+DIR+"/05_"+ sampleId +".bam", shell=True)

subprocess.call("fastqc -t "+cores+" --extract "+DIR+"/05_" + sampleId + ".bam >> "+DIR+"/fastqc_all.log 2>&1", shell=True)
# save the number of reads left in the bam file for QC plot
mappedReads = commands.getoutput("cat " + DIR + "/05_"+sampleId+"_fastqc/fastqc_data.txt | grep 'Total Sequences' | cut -f 2")
# keep files we are interested in and delete the rest
shutil.move(DIR+"/05_"+sampleId+"_fastqc.html", DIR+"/QC/fastqc_report_"+sampleId+"_from_bam.html")
shutil.move(DIR+"/05_"+sampleId+"_fastqc/Images/duplication_levels.png", DIR+"/plots/plot05_duplication_levels_from_bam_"+sampleId+".png")
os.remove(DIR+"/05_"+sampleId+"_fastqc.zip")
shutil.rmtree(DIR+"/05_"+sampleId+"_fastqc")

#####################################################################
# Write down all the read counts collected and plot the reads
# processing statistics
#####################################################################
stats01 = '\n'.join(map(str, [readsInitial, readsPostTrimming, int(mappedReads)]))
f = open(DIR+"/stats01_"+sampleId+"_reads_processing.txt", 'w')
f.write(stats01)
f.write("\n") ## adds newline at the end of the file otherwise R
              ## complains
f.close()
subprocess.call(pathto_bin+"/R_scripts/r01_read_processing_stats.R "+DIR+"/stats01_" + sampleId + 
                "_reads_processing.txt "+DIR+"/plots/plot01_read_processing_" + 
                sampleId + ".png", shell=True)

#####################################################################
# plot mapping stats (which are directly read from star final log)
#####################################################################
subprocess.call(pathto_bin+"/R_scripts/r02_mapping_stats.R "+DIR+"/star_" + sampleId + "_Log.final.out "+DIR+"/plots/plot02_mapping_" + sampleId + ".png", shell=True)

#####################################################################
# run the featureCounts for QC and make the corresponding plot
#####################################################################
# -F SAF --> specifies feature format, see http://bioinf.wehi.edu.au/featureCounts/ for details
# -O --> reads are allowed to be assigned to more than one feature (there is an explanation of why this is requested)
# -T --> number of threads
# -a --> name of the annotation file
subprocess.call(pathto_bin+"/featureCounts -F SAF -O -T 3 -a " + featureCountsDir + "/5utrs_refGene.saf -o "+
                DIR+"/5UTR_" + sampleId + "-counts.txt "+DIR+"/05_" + sampleId + ".bam >> "+
                DIR+"/logs/log06_all_featureCounts_" +sampleId+ " 2>&1", shell=True)
subprocess.call(pathto_bin+"/featureCounts -F SAF -O -T 3 -a " + featureCountsDir + "/rRNA.saf -o "+
                DIR+"/rRNA_" + sampleId + "-counts.txt "+DIR+"/05_" + sampleId + ".bam >> "+
                DIR+"/logs/log06_all_featureCounts_" +sampleId+ " 2>&1", shell=True)
subprocess.call(pathto_bin+"/featureCounts -F SAF -O -T 3 -a " + featureCountsDir + "/chrM.saf -o "+
                DIR+"/chrM_" + sampleId + "-counts.txt "+DIR+"/05_" + sampleId + ".bam >> "+
                DIR+"/logs/log06_all_featureCounts_" +sampleId+ " 2>&1", shell=True)

# plot what you just counted
subprocess.call(pathto_bin+"/R_scripts/r03_featureCounts_stats.R " + DIR + " " + sampleId + " "+
                DIR+"/plots/plot03_featureCounts_" + sampleId + ".png > "+
                DIR+"/logs/log07_featureCountsR_" + sampleId + " 2>&1", shell=True)
# since now all plots for this sample are concluded you can merge them
# in a single pdf using ImageMagick
subprocess.call("convert " + DIR + "/plots/plot0*" + sampleId + ".png "+
                DIR+"/plots/merged_plots_" + sampleId + ".pdf", shell=True)

#####################################################################
# clustering using CAGEr (remember G correction)
# temporarily frozen because CAGEr disappeared from 
# cluster in Tuebingen 
#####################################################################
## to run in Tue
subprocess.call("module load R; Rscript "+pathto_bin+"/R_scripts/r04_create_ctss_file.R " + BSgenome + 
                " "+DIR+"/05_" + sampleId + ".bam "+DIR+"/ctss_files/sample_" + sampleId + ".ctss > "+DIR+"/logs/log08_CAGEr_" + sampleId + "_ctss 2>&1", shell=True)
## ### # to run in Goe
## subprocess.call("/home/vsac/pipelines2/cage/R_scripts/r04_create_ctss_file.R " + BSgenome + " "+DIR+"/05_" + sampleId 
##                  + ".bam "+DIR+"/sample_" + sampleId + ".ctss > "+DIR+"/logs/log07_CAGEr_" 
##                  + sampleId + "_ctss 2>&1", shell=True)


#####################################################################
# clean-up (rm unnecessary files, move the rest in clear folders) ->
# this will be moved to the launcher script; everything can be deleted
# with the exclusion of ctss files and QC-related things
#####################################################################
shutil.move(DIR+"/star_" + sampleId + "_Log.final.out", DIR+"/stats")
shutil.move(DIR+"/stats01_" + sampleId + "_reads_processing.txt", DIR+"/stats")

shutil.move(DIR+"/01_" + sampleId + "-trimmed.log", DIR+"/logs/log01_" + sampleId + "-trimmed.log")
shutil.move(DIR+"/02_" + sampleId + "-trimmed.log", DIR+"/logs/log02_" + sampleId + "-trimmed.log")
shutil.move(DIR+"/03_finalTrim_" + sampleId + "-trimmed.log", DIR+"/logs/log03_finalTrim_" + sampleId + "-trimmed.log")
shutil.move(DIR + "/log_tagdust_"+ sampleId, DIR+"/logs/log04_tagdust_" + sampleId)

shutil.move(DIR+"/star_" + sampleId + "_Log.out", DIR+"/logs/log05_star_" + sampleId + "_Log.out")

shutil.move(DIR+"/05_" + sampleId + ".bam", DIR+"/bam_files/" + sampleId + ".bam")

os.remove(DIR+"/barcode_" + barcode + "_CAGE_artefacts_201507.txt")


#subprocess.call("mv "+DIR+"/log0*_" + sampleId + "* "+DIR+"/logs", shell=True)

## all that is now moved to the folder "needless_files" can be in fact removed
subprocess.call("mv "+DIR+"/*" + sampleId + "*.fastq "+DIR+"/needless_files", shell=True)
subprocess.call("mv "+DIR+"/*" + sampleId + "-counts.txt.summary "+DIR+"/needless_files", shell=True)
subprocess.call("mv "+DIR+"/*_" + sampleId + "* "+DIR+"/needless_files", shell=True)


#####################################################################
# check if all jobs are done: if not execution is interrupted for this
# sample; if yes the actions that need to be performed only once are
# performed (cleanup and stuff like this)
#####################################################################
complete = ' '.join(map(str, ["Sample", sampleId, "completed", "\n"]))
f = open(DIR+"/jobs_completed.txt", 'a')
f.write(complete)
f.close()

## the number of samples completed
check = int(commands.getoutput("wc -l "+DIR+"/jobs_completed.txt | cut -f 1 -d ' '"))
if (check < number_samples):
    check = int(commands.getoutput("wc -l "+DIR+"/jobs_completed.txt | cut -f 1 -d ' '"))
    print("The sample "+sampleId+" is NOT the last one. The pipeline just for this sample will be interrupted")
    sys.exit(0)

## when all jobs are successfully completed, perform actions that need
## to be run only once. First create the barcode splitting diagnostic
## plot
subprocess.call(pathto_bin+"/R_scripts/r00_barcode_splitting.R "+DIR+"/stats00_barcode_splitting.txt "+DIR+"/plots/plot00_barcode_splitting_barplot.png", shell=True)

shutil.move(DIR+"/fastqc_all.log", DIR+"/logs")
shutil.move(DIR+"/stats00_barcode_splitting.txt", DIR+"/stats/")

for filename in glob.glob(DIR+"/*.std.*"):
    shutil.move(filename, DIR+"/logs/")

os.remove(DIR+"/jobs_completed.txt")
os.remove(DIR+"/00_cage.fastq.gz")

## these files need to be cancelled in the "real" version
shutil.move(DIR+"/CAGE_artefacts_201507.txt", DIR+"/input/")
shutil.move(DIR+"/barcodes.txt", DIR+"/input/")
for filename in glob.glob(DIR+"/*fastq.gz"):
    shutil.move(filename, DIR+"/input/")



