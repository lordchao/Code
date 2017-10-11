source("http://bioconductor.org/biocLite.R")
biocLite("rcellminer")
biocLite("rcellminerData")
library(rcellminer)
library(rcellminerData)
searchForNscs("nib$")  
# Get Cellminer data
drugAct <- exprs(getAct(rcellminerData::drugData))
molData <- getMolDataMatrices()
# One drug
nsc <- "94600"
plots <- c("drug")
plotCellMiner(drugAct, molData, plots, nsc, NULL)
