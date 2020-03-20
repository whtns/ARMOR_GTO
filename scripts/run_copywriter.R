args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(CopywriteR)
  library(CopyhelpeR)
  library(DNAcopy)
  library(BiocParallel)
  library('tidyverse')
  library('fs')
  library('rprojroot')
  library('glue')
  library(tictoc)
})

# if (is.null(threads)) threads = 6
# if (is.null(copywriter_output_dir)) mypath = "../../output/copywriter"
# if (is.null(sample_files))   sample_files <- fs::dir_ls(mypath, glob = "*41-CL_1_recalibrated_10_intronic.bam")
# if (is.null(control_files))   control_files <- sample_files 
# if (is.null(bin_size)) bin_size = 50000
# if (is.null(suffix)) suffix = ""
# if (is.null(copywriter_capture_regions)) copywriter_capture_regions = "~/rb_pipeline/corrected_agilent_regions.bed" #must be a .bed file!

bin_size <- as.numeric(bin_size)
print(bin_size)

sample_files <- fs::dir_ls(input_dir, glob = paste0("*", samples_pattern), recurse = TRUE)
control_files <- sample_files

#BiocParallel
bp.param <- MulticoreParam(workers = 6)

## ------------------------------------------------------------------------

humanreadable <- function(bin_size){
  bin_size <- bin_size/1000
  bin_size <- paste0(toString(bin_size), "kb")
}

run_copywriter <- function(bin_size, sample_files, control_files, out_dir, ...){
  # browser()
  
  if(length(sample_files) == 0 | !is_file(sample_files)){
    error("please provide sample_files that exist")
  }
  
  dir.create(out_dir)
  
  #preCopywriteR
  
  # data.folder <- tools::file_path_as_absolute(file.path(in_dir))
  preCopywriteR(output.folder = file.path(out_dir), bin.size = bin_size, ref.genome = "hg38", "")
  

  
  if(!all(sample_files == control_files)){
    sample_files <- c(sample_files, control_files)
    
    control_files <- c(control_files, control_files)
  }
  

  sample.control <- data.frame(sample = c(sample_files), 
                               control = c(control_files))
  
  CopywriteR(
    sample.control = sample.control,
    destination.folder = file.path(out_dir),
    reference.folder = file.path(out_dir, paste0("hg38_", humanreadable(bin_size))),
    bp.param = bp.param,
    ...
  )
  
  #segment and visualize results
  plotCNA(destination.folder = file.path(out_dir))
  
}

## ------------------------------------------------------------------------

# debug(CopywriteR::plotCNA)

if (dir.exists(copywriter_output_dir)) fs::dir_delete(copywriter_output_dir)

tic()
run_copywriter(bin_size, sample_files, control_files, copywriter_output_dir, keep.intermediary.files = T)
toc()
