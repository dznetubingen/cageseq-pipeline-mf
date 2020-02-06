#!/usr/bin/Rscript --vanilla

##################################################

## upload bam file, write down a ctss file using CAGEr package. A ctss
## file is characterized by the following format (tab separated):
## chorm start strand count

##################################################

options(stringsAsFactors=FALSE)

## three arguments: genome, input bam file name, output ctss file name
args <- commandArgs(trailingOnly = TRUE)

## assign command line arguments to variables
bsGenome = args[1]
inBamFile = args[2]  ## file bam file created by the pipeline
outCtssName = args[3] ## output file to store the ctss file

## If there aren't exactly three arguments provided in command line
## display help and exit
if(length(args) != 3) {
        args <- c("--help")
}
## Help section
if ("-h" %in% args || "--help" %in% args) {
        cat("
        Please pass 3 arguments to the script:
	BS Genome name: to be used by CAGEr for mapping
        input file (bam): name of the input bam file
        output name: name of the ctss file created as output

        Example:
        r04_create_ctss_file.R BSgenome.Hsapiens.UCSC.hg19 10196_caudate.bam sample_10196_caudate.ctss\n\n")
        q(save="no")
}

############################
## start actual execution
############################

## load required libraries
library(CAGEr)
library(bsGenome, character.only=TRUE)

## define CAGEset object
myCAGEset <- new("CAGEset",
                 genomeName = bsGenome,
                 inputFiles = inBamFile,
                 inputFilesType = "bam",
                 sampleLabels = "sample")

## load the CAGE data
## updated 2015/09/07 to obtain "full G correction"
## NB: default behaviour is to exclude reads with mapping quality < 20
getCTSS(myCAGEset, removeFirstG = TRUE, correctSystematicG = TRUE)

## create ctss file
ctss <- CTSStagCount(myCAGEset)

## force the coordinate column to be integer, otherwise it creates
## problems in the step that creates the count table to be fed to DE
ctss[,2] <- as.integer(ctss[,2])

## save the ctss to a file with appropriate name
write.table(ctss, file=outCtssName, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

