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

  mass_tol <- 1e-4

  if(!('file'  %in% names(data))){
    data <- data %>% mutate(file="file")
  }
  if(is.na(group))
    group='group'
  if(!(group  %in% names(data))){
    data <- data %>% mutate(group="group")
  }else if(group!='group'){
    data <- data %>% mutate_(group=sprintf("`%s`", group))
  }
  message("Files to process: ")
  print(unique(data$file))
  message("Groups to process: ")
  print(unique(data$group))
  spectra <- unique(data$seqNum)
  mono_spectra <- unique((data %>% filter(peak=='0'))$seqNum)
  nomono_spectra <- setdiff(spectra, mono_spectra)
  require(multidplyr)
  cl <- create_cluster()
  cluster_library(cl, "dplyr")
  set_default_cluster(cl)

  data %>%
    filter(abs(masserror)<mass_tol) %>%
#    group_by(group,file, seqNum, rt) %>%
    partition(group,file, seqNum, rt) %>%
    do({
      dd <- .
      res <- data.frame()
      elements <- setdiff(unique(dd$peak), '0')
      for(el in elements){
        good_aa <- (dd %>% filter(peak==el & I >0) %>% distinct(ion))$ion
        i0s <- (dd %>% filter(ion %in% good_aa & peak=='0'))$I
        i1s <- (dd %>% filter(ion %in% good_aa & peak==el))$I
        ns <- (dd %>% filter(ion %in% good_aa & peak==el))$n

        rs <- i1s/i0s
        i0s <- i0s[rs<1]
        i1s <- i1s[rs<1]
        ns <- ns[rs<1]
        rs <- rs[rs<1]
        rs_good <- (abs(rs-median(rs))<1.5*mad(rs))
        r <- sum(i1s[rs_good]/ns[rs_good])/sum(i0s[rs_good])
        if(!is.na(r))
          res <- res %>% bind_rows(data.frame(peak=el, gamma=r, logI = log2(sum(i1s))))
      }
      if(nrow(res)>0){
        # res$seqNum = dd$seqNum[[1]]
        # res$rt = dd$rt[[1]]
        res$ion = 'any'
      }
      res
    }) %>% collect() %>%  ungroup() -> any_aa
  data %>%
    filter(abs(masserror)<mass_tol) %>%
    filter(n>0 & I>0) %>%
    mutate(gamma=isoratio/n, logI=log2(I)) %>%
    select(group,file,peak, gamma, logI, seqNum, rt, ion) %>%
    bind_rows(any_aa) %>% ungroup() -> all_aa

  for(pp in unique(all_aa$peak)){
    vals <- (all_aa %>% filter(peak==pp))$gamma
    vmin <- quantile(vals, 0.05)
    vmax <- quantile(vals, 0.95)
    all_aa %>%
      filter(peak==pp) %>%
      group_by(ion, file ) %>%
      summarize(m = median(gamma), md = mad(gamma)) -> mmad
    aa_data <- all_aa %>%
      filter(peak==pp) %>%
      left_join(mmad %>% ungroup(), by=c('ion','file')) %>%
      filter(md>0) %>%
      filter(gamma>(m-3*md) & gamma < (m+3*md))
    aa_data$ion <- factor(aa_data$ion, levels = (mmad %>% summarise(md=median(md)) %>% arrange(md))$ion)
    aa_data %>%
      group_by(ion, group, file) %>%
      summarize(mean=mean(gamma)*100, median=median(gamma)*100, wmean = weighted.mean(gamma, logI)*100, sd = sd(gamma)*100, se = sd(gamma)/sqrt(n())*100, n=n()) %>%
      mutate(CV = se/wmean*100) %>%
      arrange(ion) -> res_f
    res_f %>% write_csv(file.path(resultPath, sprintf("file_summary_%s.csv", pp)))
    res_f %>%
      summarize(gmean=mean(mean), gmedian=mean(median), gwmean=mean(wmean), sd_mean=sd(mean), sd_median=sd(median), sd_wmean=sd(wmean), n=n()) %>%
      mutate(
        cv_mean = sd_mean/gmean/sqrt(n)*100,
        cv_median = sd_median/gmedian/sqrt(n)*100,
        cv_wmean = sd_wmean/gwmean/sqrt(n)*100
      ) %>%
      arrange(cv_wmean) %>%
      write_csv(file.path(resultPath, sprintf("group_summary_%s.csv", pp)))
  }
  if(length(nomono_spectra)>0){
    cN <- (data %>% ungroup() %>% filter(seqNum %in% nomono_spectra) %>% filter(peak=='13C') %>% slice(1:1))$n
    gC <- 0.01
    data %>%
      filter(seqNum %in% nomono_spectra) %>%
      # filter(abs(masserror)<mass_tol) %>%
      filter(n>0 & I>0) %>%
      mutate(gamma=isoratio/n*cN*gC, logI=log2(I)) -> narrow_data
    narrow_data %>%
      group_by(peak, file ) %>%
      summarize(m = median(gamma), md = mad(gamma)) -> narrowmmad

    narrow_data%>%
      left_join(narrowmmad %>% ungroup(), by=c('peak','file')) %>%
      filter(md>0) %>%
      filter(gamma>(m-3*md) & gamma < (m+3*md)) %>%
      group_by(ion,peak, group, file) %>%
      summarize(mean=mean(gamma)*100, median=median(gamma)*100, wmean = weighted.mean(gamma, logI)*100, sd = sd(gamma)*100, se = sd(gamma)/sqrt(n())*100, n=n()) %>%
      mutate(CV = se/wmean*100) %>%
      arrange(ion) -> narrowres_f
    narrowres_f %>% write_csv(file.path(resultPath, sprintf("narrow_file_summary.csv")))
    narrowres_f %>%
      summarize(gmean=mean(mean), gmedian=mean(median), gwmean=mean(wmean), sd_mean=sd(mean), sd_median=sd(median), sd_wmean=sd(wmean), n=n()) %>%
      mutate(
        cv_mean = sd_mean/gmean/sqrt(n)*100,
        cv_median = sd_median/gmedian/sqrt(n)*100,
        cv_wmean = sd_wmean/gwmean/sqrt(n)*100
      ) %>%
      arrange(cv_wmean) %>%
      write_csv(file.path(resultPath, sprintf("narrow_group_summary.csv")))

  }

  gg <- all_aa %>% distinct(group)
  rmarkdown::render(system.file("Rmd/isoMS_report.Rmd", package =getPackageName()), envir= sys.frame(sys.nframe()), output_file=file.path(resultPath,"isoMS_report.html"))
}
