#!/usr/bin/Rscript --vanilla

options(stringsAsFactors=FALSE)

## parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

## assign command line arguments to variables
inFile = args[1]  ## file "stats00_barcode_splitting.txt"
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
        inFile: stats00_barcode_splitting.txt
        output name: name of the png created as output

        Example:
        r00_barcode_splitting.R stats00_barcode_splitting.txt plot00_barcode_splitting_barplot.png\n\n")
        q(save="no")
}

############################
## start actual execution
############################
# read table containg barcode splitting counts
d <- read.table(args[1])
l <- length(d[,1])

# counts are stored in second column; I skip the first row because it
# contains the total number of reads
counts <- d[2:l,2]
# what each number represents is written in the first column
names(counts) <- d[2:l,1]

# the last "entry" represents the number of unmatched reads, given by the total number
# of reads minus all the reads assigned to a given barcode
counts[l] <- d[1,2] - sum(counts) 
names(counts)[l] <- "unmatched"

# calculate percentages, which is what I need for the plot
perc <- counts/sum(counts)*100
ylim <- c(0, 1.1*max(perc))

# store plot in a png file
png(args[2])
par(las=2, mar=c(7, 4, 4, 2))

# Plot, and store x-coordinates of bars in xx to write the number of reads in the 
# right place of the graph
xx <- barplot(perc, width = 0.85,
              main=paste("Barcode splitting statistics\nTotal reads ", d[1,2], sep=""),
              ylim=ylim, ylab="Percentage (%)")
# Add actual number of reads per bar at top of bars
text(x = xx, y = perc, label = counts, pos = 3, cex = 0.8, col = "red")

# this avoids getting a couple of lines of message by running the script with Rscript
garbage <- dev.off()


