processExperimentDesign <- function(){
  library(isoms)
  args_ <- commandArgs(trailingOnly = TRUE)
  expDesignFile <- args_[[1]]
  if(!file.exists(expDesignFile)){
    message("You must provide corrent experimentDesign file location")
    return(0)
  }
  expDesign <- suppressMessages(read_csv(expDesignFile))
  message("Samples to analyze: ", toString(unique(expDesign$Sample)))
  message("Number of controls: ", nrow(expDesign %>% filter(Sample=='Control')))

  expDesign %>%
    rowwise() %>%
    do({
      xx = .
      fd <- suppressWarnings(suppressMessages(read_csv(xx$File)))
      if("totIonCurrent" %in% names(fd))
        fd <- fd %>% mutate(tic=totIonCurrent) %>% select(-totIonCurrent)
      fd %>%
      mutate(file=xx$File, group=xx$Sample, loading=xx$Loading)
    }) -> data_
  summarizeImmoniums(
    data = data_,
    group='group',
    resultPath = dirname(expDesignFile)
  )
}

# ff <- list.files("~/tmp/Immonium/Peptome/20170901/", pattern = ".*_\\d+.*_gaussian2.csv")
# data <- bind_rows(
#   lapply(
#     ff,
#     function(x){
#       print(x)
#       read_csv(file.path("~/tmp/Immonium/Peptome/20170901", x)) %>%
#         mutate(
#           file=sub(".*Imm_(.+)_gaussian.*","\\1",x),
#           group=sub("_\\d+","",file),
#           loading=as.numeric(sub(".*_(\\d+)_gaussian.*","\\1", x)))
#     }
#   )
# )
#
#
