#!/usr/bin/Rscript --vanilla

## load needed libraries and set options
library(ggplot2)
options(stringsAsFactors=FALSE)

## parse command line arguments
args <- commandArgs(TRUE)

## assign command line arguments to variables
inFileName = args[1]  ## file "stats01_03121_front_reads_processing.txt"
outPNGName = args[2] ## output file name for plotting

## If there aren't exactly two arguments provided in command line display help
## and exit
if(length(args) != 2) {
        args <- c("--help")
}
## Help section
if ("-h" %in% args || "--help" %in% args) {
        cat("
        Please pass 2 arguments to the script:
        in file name: stats01_<sampleId>_reads_processing.txt
        output name: name of the png created as output

        Example:
        r01_read_processing_stats.R stats01_03121_front_reads_processing.txt plot01_read_processing_03121_front.png\n\n")
        q(save="no")
}

############################
## start actual execution
############################

# read the file containing mapping stats
d <- read.table(inFileName)
# from it calculate the percentages
percentages <- d[,1]/d[1,1] *100
## create dataframe for ggplot2 plotting
processingStage = c("Starting reads", "After trimming", "Uniquely mapping")
data <- data.frame(processingStage=factor(processingStage, levels=processingStage), 
                   Percentage = percentages)
## plot (on the png which name is provided as input)
png(outPNGName)
ggplot(data = data, aes(x = processingStage, y=Percentage, fill=processingStage)) +  
        geom_bar(stat="identity") + 
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
        labs(title="Retained reads")
## this avoids getting a couple of lines of message by running the script with Rscript
garbage <- dev.off()