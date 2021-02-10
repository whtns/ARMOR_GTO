args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library('tidyverse')
  library('fs')
  library('rprojroot')
  library(ichorCNA)
})

print(outDir)
print(samplewig)

ichorCNA::runIchorCNA(samplewig, 
            gcWig = fs::path_package("ichorCNA", "extdata/gc_hg19_1000kb.wig"), 
            mapWig = fs::path_package("ichorCNA", "extdata/map_hg19_1000kb.wig"))

sessionInfo()
date()


