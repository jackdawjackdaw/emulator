#
# plot the model points and 
# the splits too

#library(Hmisc) # for errbar
library(gtools)

cols <- c("black", "coral", "blue", "green", "darkorchid", "cornflowerblue", "orange",  "chocolate1", "cadetblue", "red")

#modelData <- read.table("model-cut.dat") # temp, mean, var
modelData <- read.table("breaksit.dat") # temp, mean, var

plot(modelData[,1], modelData[,2], col=cols[1], xlab="temp (arb)", ylab="mag (arb)", ylim=range(0,1.2))

# read the output folder
inFiles <- list.files("output/",full.names=TRUE)
inFiles <- rev(mixedsort(inFiles)) # sort the right way around

for(i in 1:length(inFiles)){
  splitData <- read.table(inFiles[i])
  lines(splitData[,1], splitData[,2], col=cols[i+1])
}
               
