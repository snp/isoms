mzMLtoCSV <- function(pattern="*.mzML"){
  library(isoms)
  args_ <- commandArgs(trailingOnly = TRUE)
  dir_ = "."
  if(length(args_) > 0)
    dir_ <- args_[[1]]
  if(length(args_) > 1)
    pattern <- args_[[2]]
  if(!dir.exists(dir_)){
    message("First argument has to be the directory containing mzML files")
    return(0)
  }
  ff <- list.files(path=dir_, pattern = pattern)
  if(length(ff)<1){
    message("Now mzML files found in directory [", dir_,"] using pattern [", pattern, "].")
    return(0)
  }
  for(f in ff){
    in_f <- file.path(dir_, f)
    outf <- sub(".mzML$","_fit.csv", in_f)
    if(!file.exists(outf)){
      message("Converting file [", f, "]")
      analyze_immoniums(file = in_f, width = 0.0015) %>%
        write_csv(outf)
    }else{
      message("File already processed: [", f, "].")
    }
  }
}
