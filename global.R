# Author: Matt Watts
# Date: 12 Dec 2014
# Purpose: RunMarxanLPILPTas web app global.R

library(shiny)
library(PBSmapping)
library(maptools)
library(sp)

cat(paste0("hello\n"))
cat(paste0(getwd(),"\n"))

sMarxanDir <- getwd()

# find how many runs from input.dat
inputdat <- readLines(paste(sMarxanDir,"/input.dat",sep=""))
iParam <- which(regexpr("NUMREPS",inputdat)==1)
iNUMREPS <<- as.integer(unlist(strsplit(inputdat[iParam], split=" "))[2])

irefreshinput <<- 0
isavetargetspf <<- 0
