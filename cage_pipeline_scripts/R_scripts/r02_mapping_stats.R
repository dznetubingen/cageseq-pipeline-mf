#!/usr/bin/Rscript --vanilla

library(ggplot2)
options(stringsAsFactors=FALSE)

## parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
## assign command line arguments to variables
inFile = args[1]  ## file "..Log.final.out" that STAR outputs
outPNGName = args[2] ## output file for plotting
## If there aren't exactly two arguments provided in command line display help
## and exit
if(length(args) != 2) {
        args <- c("--help")
}
## Help section
if ("-h" %in% args || "--help" %in% args) {
        cat("
        Please pass 2 arguments to the script:
        inFile: file something.Log.final.out that STAR outputs
        output name: name of the png created as output

        Example:
        r02_mapping_stats.R star_10196_caudate_Log.final.out mapping_output_plot.png\n\n")
        q(save="no")
}

############################
## start actual execution
############################

## read input file
d <- read.table(inFile, fill=TRUE, sep="\t")
## extract numbers we are interest in, parse where needed to remove unwanted characters
perc <- as.vector(d[,2])[c(9,26,28:30)]
perc <- as.numeric(gsub("%", "", perc))
input_reads <- as.numeric(as.vector(d[,2])[5])
## labels, needed to create the legend of the barplot
perc_labels <- c("Unique mappers", "Mapping > 1 locus", 
                 "Unmapped: too many\nmismatches", "Unmapped: too short",
                 "Unmapped: other")
## create dataframe necessary for ggplot2 plotting: first column is the mapping
## class, second column is the actual percentage for the corresponding mapping
## class
data <- data.frame(mappingClass = factor(perc_labels, levels = perc_labels), 
                   Percentage = perc)
## actual plotting using ggsave functionality
ggsave(outPNGName, dpi=300, width = 13, height = 13, units = "cm",
       ggplot(data = data, aes(x = mappingClass, y=Percentage, fill=mappingClass)) + 
               geom_bar(stat="identity") + 
               geom_abline(slope = 0, intercept = 80, lty= "dashed") + 
               theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
               labs(title=paste("Mapping statistics\nInput reads ", input_reads, sep="")) +
               ylim(c(0,100))
)
