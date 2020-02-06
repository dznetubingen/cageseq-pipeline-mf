#!/usr/bin/Rscript --vanilla


library(ggplot2)
options(stringsAsFactors=FALSE)

## parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

## we require exactly three arguments (args[1] is the directory where
## we are running, args[2] is the sampleId, used to make input file
## names; args[3] is the output png file for plotting).

## assign command line arguments to variables
dir =  args[1]
sampleId = args[2]
outPNGName = args[3]
## If there aren't two arguments display help and exit
if(length(args) != 3) {
        args <- c("--help")
}
## Help section
if ("-h" %in% args || "--help" %in% args) {
cat("
        Please pass 3 arguments to the script:
	directory: where counts summaries are to be found
        sample Id: sample ID as given in barcode file;
        output name: name of the png created as output

        Example:
        r03_featureCounts_stats.R /home/user/run1/ 10196_caudate output_plot.png\n\n")
        q(save="no")
}

############################
## start actual execution
############################

## types of featureCounts inputs I need to read
features <- c("5UTR","rRNA","chrM")
suffix <- c("-counts.txt.summary")
## create the names for all input file to be read
inputNames <- paste(features, sampleId, sep="_")
inputNames <- paste(dir, inputNames, suffix, sep="")
## create vector to store the percentage of reads mapping to each of the features
fractions <- c()
## open all files and read the numbers of interest
for (i in 1:length(inputNames)){
        d <- read.table(inputNames[i], header=TRUE)
        assigned <- d[1,2]
	unassigned <- d[4,2]
        total <- assigned + unassigned
        fractions <- c(fractions, assigned/total)
}
## create dataframe for ggplot2 plotting
data <- data.frame(Features=factor(features, levels=features), 
                   Percentage = fractions*100)
## plot (on the png which name is provided as input)
ggsave(outPNGName, dpi=300, width = 11, height = 11, units = "cm",
       ggplot(data = data, aes(x = Features, y=Percentage, fill=Features)) +  
               geom_bar(stat="identity") + 
               theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
               geom_abline(slope = 0, intercept = 60, lty= "dashed") +
               geom_abline(slope = 0, intercept = 10, lty= "dashed") +
               geom_abline(slope = 0, intercept = 5, lty= "dashed") +
               ylim(c(0,100))  + 
               labs(title="Tag distribution across features")
)
