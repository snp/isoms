#' summarizeImmoniums - generate the summary of the immonium analysis
#' Given the results of \code{analyzeImmoniums} function, produces a nice html
#' report
#'
#' @param data - a dataframe returned by \code{analyzeImmoniums} function.
#' If \code{NA} then using \code{files}
#' @param files - vector of csv-files containing results of \code{analyzeImmoniums}
#' @param group - name of the column in data
#'
#' @return
#' @export
#'
#' @examples
summarizeImmoniums <- function(data=NA, files=NA, group=ifelse(is.na(data), 'file', NA), resultPath="./isoMS_result"){
  if(is.na(data)){
    if(class(files) == "character")
      data <- bind_rows(
        lapply(files, function(f){
          if(file.exists(f))
            read_csv(f) %>% mutate(file=f)
          else
            data.frame()
        })
      )
  }
  if(nrow(data)<2){
    warning("You have to provide either data frame or vecor of file names to proceed")
    return(0)
  }
  if (!dir.exists(resultPath))
    dir.create(resultPath, recursive = TRUE)
  render(system.file("Rmd/isoMS_report.Rmd", package =getPackageName()), envir= sys.frame(sys.nframe()), output_file=file.path(resultPath,"isoMS_report.html"))
}
