#!/usr/bin/env Rscript
## script to create consensus clusters across distinct 

## r05_make_count_table_CAGEr BSgenome.Hsapiens.UCSC.hg38 input_dir

## requires arguments:
## arg[1] -> BSgenome corresponding to organism
## arg[2] -> input_data_directory; where the ctss files are stored

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
  args <- c("--help")
}

## Help section
if("-h" %in% args ||"--help" %in% args) {
  cat("
        Please pass 2 arguments to the script:
        BSgenome db: database;
	path_to_input: directory where you will store the input files;
 
        Example:
        make_count_table_CAGEr.R BSgenome.Hsapiens.UCSC.hg19 input_dir\n\n")
        q(save="no")
}


reference=args[1]
input_dir <- args[2]

######################################
#reference="BSgenome.Hsapiens.UCSC.hg19"
#input_dir <- "ctss_files"
######################################

options(stringsAsFactors=FALSE)
library(CAGEr)
library(reference, character.only=TRUE)

## load all ctss files in the specified input directory
pathsToInputs <- list.files(input_dir, full.names=TRUE)
sample_labels <- list.files(input_dir)
sample_labels <- gsub(".ctss", "", sample_labels)

print("Set cageset object")
myCAGEset <- new("CAGEset", genomeName = reference,
                 inputFiles = pathsToInputs,
                 inputFilesType = "ctss", sampleLabels = sample_labels)
getCTSS(myCAGEset, removeFirstG = FALSE)

print("Setting ctss object")
ctss <- CTSStagCount(myCAGEset)

print("Normalize counts")
normalizeTagCount(myCAGEset, method="none")

print("clustering...")
clusterCTSS(object = myCAGEset, threshold = 10, thresholdIsTpm = FALSE,
         nrPassThreshold = 2, method = "distclu", maxDist =20,
         removeSingletons = TRUE)

print("aggreagate cluters")
aggregateTagClusters(myCAGEset, tpmThreshold = 0, qLow = NULL, qUp = NULL, maxDist = 100)

countTable <- myCAGEset@consensusClustersTpmMatrix
head(countTable)

for_deseq_consensus_cluster <- consensusClusters(myCAGEset)

for_deseq_consensus_cluster$id <- paste(for_deseq_consensus_cluster$chr,
                               for_deseq_consensus_cluster$start, for_deseq_consensus_cluster$end,
                               for_deseq_consensus_cluster$strand, sep="_")

rownames(countTable) <- for_deseq_consensus_cluster$id

write.table(countTable, file=paste(input_dir,"/","countTable.txt", sep=""), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
#write.table(colnames(countTable), file="sampleIDs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)
