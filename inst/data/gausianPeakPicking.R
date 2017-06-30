#source("https://bioconductor.org/biocLite.R")
#biocLite("mzR")
#library(RforProteomics)
library(mzR)
library(parallel)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(readr)

elements <- data_frame(
  element = c("C","N","O","H"),
  mass = c(12., 14.00307400443, 15.99491461957,1.00782503223),
  delta = c(
    13.00335483507 - 12.0,
    15.00010889888 - 14.00307400443,
    17.999161 - 15.99491461957,
    2.01410177812 - 1.00782503223
  ),
  gamma = c(
    0.0107 / 0.9893,
    0.00364 / 0.99636,
    0.00205 / 0.99757,
    0.000115/0.999885
  )
)

immoniumIons = list(
  'P' = c('C' = 4, 'H' = 8, 'N' = 1), # 70.06587
  'L' = c('C' = 5, 'H' = 12, 'N' = 1), # 86.09707
  'R' = c('C' = 5, 'H' = 13, 'N' = 4),  # 129.11358
  'K' = c('C' = 5, 'H' = 13, 'N' = 2),  # 101.10784
  'K(CH2O)' = c('C' = 6, 'H' = 15, 'N' = 2, 'O' = 2),
  'A' = c('C'= 2, 'H' = 6, 'N' = 1),   # 46...
  'H' = c('C'= 5, 'H' = 8, 'N' = 3),   # 110.07174
  'F' = c('C'= 8, 'H' = 10, 'N'= 1),   # 120.08115
  'S' = c('C'= 2, 'H' = 6, 'N'= 1, 'O' = 1),  # 60.04514
  'V' = c('C' = 4, 'H' = 10, 'N' = 1), # 72.08151
  'Y' = c('C' = 8, 'H' = 10, 'N' = 1, 'O' = 1),    # 136.07593
  'N' = c('C' = 3, 'H' = 7, 'N' = 2, 'O' = 1),
  'D' = c('C' = 3, 'H' = 6, 'N' = 1, 'O' = 2),
  'Q' = c('C' = 4, 'H' = 9, 'N' = 2, 'O' = 1),
  'E' = c('C' = 4, 'H' = 8, 'N' = 1, 'O' = 2),
  'T' = c('C' = 3, 'H' = 8, 'N' = 1, 'O' = 1),
  'Pyr' = c('C' = 4, 'H' = 6, 'N' = 1, 'O' = 1))

monoMass <- function(comp, charge=0){
  elements <- data_frame(
    element = c("C","N","O","H"),
    mass = c(12., 14.00307400443, 15.99491461957,1.00782503223),
    delta = c(
      13.00335483507 - 12.0,
      15.00010889888 - 14.00307400443,
      17.999161 - 15.99491461957,
      2.01410177812 - 1.00782503223
    ),
    gamma = c(
      0.0107 / 0.9893,
      0.00364 / 0.99636,
      0.00205 / 0.99757,
      0.000115/0.999885
    )
  )
  mass <- as.numeric( charge * elements$mass[elements$element=='H'])
  for(el in names(comp)){
    if(el %in% elements$element)
      mass <- mass + comp[el] * elements$mass[elements$element==el]
    else
      warning(sprintf('Unknown element: %s', el))
  }
  unname(mass)
}

rSquared <- function(model, y) {
  ## Determine R2 of a fitted model
  if (class(model)!="try-error"){
    ssTot <- sum((y - mean(y, na.rm=TRUE))^2, na.rm=TRUE)
    ssRes <- sum( (y - predict(model))^2 , na.rm=TRUE)
    r2 <- 1 - ssRes/ssTot
  } else {
    r2 <- NA_real_
  }
  return(r2)
}

fitPeak <- function(x, y, startPars=c('mu' = 70.065, 'I' = 5e3, 'sigma' = 3e-4),maxAttempts=1000){
  fctStr <- "I/sigma*exp(-1/2*(x-mu)^2/sigma^2)"
  fitFct <- as.formula(paste("y ~", fctStr))
  xVec <- x
  yVec <- y
  varyPars <- 0
  attempts <- 0
  repeatLoop <- TRUE
  ## Check if number of non-missing values is sufficient
  ## (NLS can only handle data with at least three non-missing values)
  validValues <- !is.na(yVec)
  if (sum(validValues) <=2){
    m <- NA
    class(m) <- "try-error"
  } else{
    ## Perform fit
    while(repeatLoop & attempts < maxAttempts){
      parTmp <- startPars * (1 + varyPars*(runif(1, -0.1, 0.1)))
      m <- try(nls(formula=fitFct, start=parTmp, data=list(x=xVec, y=yVec), na.action=na.exclude,
                   algorithm="port"),
               silent=TRUE)
      attempts <- attempts + 1
      varyPars <- 1
      if (class(m)!="try-error") {
        mcoef <- try(tidy(m))
        if(class(mcoef)!="try-error")
          repeatLoop <- FALSE
      }
    }
  }
  return(m)
}

get_isopeaks <- function(ss, mz, width=0.002, npoint=10){
  elements <- data_frame(
    element = c("C","N","O","H"),
    mass = c(12., 14.00307400443, 15.99491461957,1.00782503223),
    delta = c(
      13.00335483507 - 12.0,
      15.00010889888 - 14.00307400443,
      17.999161 - 15.99491461957,
      2.01410177812 - 1.00782503223
    ),
    gamma = c(
      0.0107 / 0.9893,
      0.00364 / 0.99636,
      0.00205 / 0.99757,
      0.000115/0.999885
    )
  )

  # Checking input data
  if(length(dim(ss)) != 2)
    return(data.frame())
  # here we'll put our result
  res <- data.frame()
  # Fitting monoisotopic peak. We allow wide window for search: 3 times width
  # finding mz-interval
  p_int <- findInterval(c(mz-3*width, mz+3*width), ss[,1])
  if(diff(p_int) < npoint)
    return(data.frame())
  sub_ss <- ss[p_int[1]:p_int[2],]
  istart = iend = which.max(sub_ss[,2])
  while(sub_ss[istart,2]>0 & istart > 1)
    istart <- istart - 1
  while(sub_ss[iend,2]>0 & iend < nrow(sub_ss))
    iend <- iend + 1
  if((iend-istart)<npoint)
    return(data.frame())
  sub_ss <- sub_ss[istart:iend,]
  # adjusting initial parameters for gaussian fit
  startPars <- c(
    "mu" = sub_ss[which.max(sub_ss[,2]),1],
    "sigma" = sd(sub_ss[sub_ss[,2]>0, 1])/2,
    "I" = sub_ss[which.max(sub_ss[,2]),2] * sd(sub_ss[sub_ss[,2]>0, 1]) / 2
  )
  model = fitPeak(sub_ss[,1], sub_ss[,2], startPars = startPars, maxAttempts = 100)
  r2 <- rSquared(model, sub_ss[,2])
  # plot(sub_ss)
  # xx = seq(min(sub_ss[,1]), max(sub_ss[,1]),length.out=400)
  # yy = predict(model, newdata=list(x=xx))
  # lines(xx,yy)
  if(class(model)!='try-error' & r2 > 0.9){
    v <- coef(model)
    res <- tidy(model) %>% mutate(peak='0')
    # 13C
    p_int <- findInterval(
      c(v['mu'] + elements$delta[elements$element=='C']-width,
        v['mu'] + elements$delta[elements$element=='C']+width),
      ss[,1])
    if(diff(p_int)>npoint){
      sub_ss <- ss[p_int[1]:p_int[2],]
      istart = iend = which.max(sub_ss[,2])
      while(sub_ss[istart,2]>0 & istart > 1)
        istart <- istart - 1
      while(sub_ss[iend,2]>0 & iend < nrow(sub_ss))
        iend <- iend + 1
      if((iend-istart+1) > npoint){
        sub_ss <- sub_ss[istart:iend,]
        startPars <- c(
          v['mu'] + elements$delta[elements$element=='C'],
          v['sigma'],
          v['I'] * elements$gamma[elements$element=='C'] * 4)
        model_C = fitPeak(sub_ss[,1], sub_ss[,2], startPars = startPars, maxAttempts = 100)
        r2_C = rSquared(model_C, sub_ss[,2])
        # plot(sub_ss)
        if(class(model_C)!='try-error' & r2_C > 0.5){
          # xx = seq(min(sub_ss[,1]), max(sub_ss[,1]),length.out=400)
          # yy = predict(model_C, newdata=list(x=xx))
          # lines(xx,yy)
          vC <- coef(model_C)
          res <- res %>% bind_rows(tidy(model_C) %>% mutate(peak='C'))
        }
      }
    }
    # 15N
    p_int <- findInterval(
      c(v['mu'] + elements$delta[elements$element == 'N'] - width,
        v['mu'] + elements$delta[elements$element == 'N'] + width),
      ss[,1])
    if(diff(p_int) > npoint){
      sub_ss <- ss[p_int[1]:p_int[2],]
      istart = iend = which.max(sub_ss[,2])
      while(sub_ss[istart,2]>0 & istart > 1)
        istart <- istart - 1
      while(sub_ss[iend,2]>0 & iend < nrow(sub_ss))
        iend <- iend + 1
      if((iend-istart+1) > npoint){
        sub_ss <- sub_ss[istart:iend,]
        startPars <- c(
          v['mu'] + elements$delta[elements$element=='N'],
          v['sigma'],
          v['I'] * elements$gamma[elements$element=='N'] * 4)
        model_N = fitPeak(sub_ss[,1], sub_ss[,2], startPars = startPars, maxAttempts = 100)
        r2_N <- rSquared(model_N, sub_ss[,2])
        # plot(sub_ss)
        if(class(model_N)!='try-error' & r2_N>0.5){
          # xx = seq(min(sub_ss[,1]), max(sub_ss[,1]),length.out=400)
          # yy = predict(model_N, newdata=list(x=xx))
          # lines(xx,yy)
          vN <- coef(model_N)
          res <- res %>% bind_rows(tidy(model_N) %>% mutate(peak='N'))
        }
      }
    }
    # 2H
    p_int <- findInterval(
      c(v['mu'] + elements$delta[elements$element=='H']-width,
        v['mu'] + elements$delta[elements$element=='H']+width),
      ss[,1])
    if(diff(p_int)>npoint){
      sub_ss <- ss[p_int[1]:p_int[2],]
      istart = iend = which.max(sub_ss[,2])
      while(sub_ss[istart,2]>0 & istart > 1)
        istart <- istart - 1
      while(sub_ss[iend,2]>0 & iend < nrow(sub_ss))
        iend <- iend + 1
      if((iend-istart+1) > npoint){
        sub_ss <- sub_ss[istart:iend,]

        startPars <- c(
          v['mu'] + elements$delta[elements$element=='H'],
          v['sigma'],
          v['I'] * elements$gamma[elements$element=='H'] * 4)
        model_H = fitPeak(sub_ss[,1], sub_ss[,2], startPars = startPars, maxAttempts = 100)
        r2_H <- rSquared(model_H, sub_ss[,2])
        if(class(model_H)!='try-error' & r2_H>0.5){
          vH <- coef(model_H)
          # message(vH['mu'] - v['mu'], " vs ",elements$delta[elements$element=='H'])
          res <- res %>% bind_rows(tidy(model_H) %>% mutate(peak='H'))
        }
      }
    }
    # 18O
    p_int <- findInterval(
      c(v['mu'] + elements$delta[elements$element=='O']-width,
        v['mu'] + elements$delta[elements$element=='O']+width),
      ss[,1])
    if(diff(p_int)>npoint){
      sub_ss <- ss[p_int[1]:p_int[2],]
      istart = iend = which.max(sub_ss[,2])
      while(sub_ss[istart,2]>0 & istart > 1)
        istart <- istart - 1
      while(sub_ss[iend,2]>0 & iend < nrow(sub_ss))
        iend <- iend + 1
      if((iend-istart+1) > npoint){
        startPars <- c(
              v['mu'] + elements$delta[elements$element=='O'],
              v['sigma'],
              v['I'] * elements$gamma[elements$element=='O'] * 4)
        model_O = fitPeak(sub_ss[,1], sub_ss[,2], startPars = startPars, maxAttempts = 100)
        r2_O <- rSquared(model_O, sub_ss[,2])
        if(class(model_O)!='try-error' & r2_O>0.5){
          vO <- coef(model_O)
          # message(vH['mu'] - v['mu'], " vs ",elements$delta[elements$element=='H'])
          res <- res %>% bind_rows(tidy(model_O) %>% mutate(peak='O'))
        }
      }
    }
  }
  return (res)
}

get_immoniums <- function(ss, ions = c('P','L','V'), width=0.002, npoint=5){
  result <- data.frame()
  for(ion in ions){
    if(!(ion %in% names(immoniumIons))){
      warning(sprintf("Don't know ion [%s]", ion))
      return(NA)
    }
    mz <- monoMass(immoniumIons[[ion]])
    isopeaks <- get_isopeaks(ss, mz, width=width, npoint = npoint)
    if(nrow(isopeaks)>3)
      result <- result %>%
        bind_rows(isopeaks %>% mutate(aa=ion))
  }
  result
}

analyze_immoniums <- function(file='/Volumes/Seagate 2T/Immonium/Fusion/HeLa/Conc/20170408_F_04_HeLa_Imm_200ng.mzML'){
  message(sprintf("Reading file [%s]", file))
  msrun <- openMSfile(file)
  hd <- header(msrun)
  immscans <- which(hd$collisionEnergy>49 & (hd$msLevel>1) & (hd$lowMZ<71))
  message(sprintf("\t%d immonium scans found", length(immscans)))
  result <- data.frame()
  message("Processing:")
  pb <- txtProgressBar(min=0, max=length(immscans))
  idx = 0

  result <- result %>%
    bind_rows(do.call(rbind, lapply(
      immscans, function(ii){
      ss <- peaks(msrun, ii)
      hh <- header(msrun, ii)
      idx <<- idx + 1
      setTxtProgressBar(pb, idx)
      get_immoniums(ss, c('L','P','V','N','T','H','Pyr'), width = 0.0016) %>%
        filter(is.numeric(estimate)) %>%
        mutate(
          rt = hh$retentionTime,
          scan = hh$seqNum)
    })#, mc.cores=ncpu)
  ))

  result$n <- apply(result, 1, function(x)immoniumIons[[x["aa"]]][x["peak"]])

  result %>%
    select(term, estimate, peak, aa, rt, scan, n) %>%
    group_by(scan, rt) %>%
    do({
      dd <- .
      res <- data.frame()
      elements <- setdiff(unique(dd$peak), '0')
      el_aa <- list(
        'C' = c('L','P','V'),
        'N' = c('L','P','V'),
        'H' = c('P','V', 'L'),
        'O' = c('Pyr','T'))
      for(el in elements){
        good_aa <- (dd %>% filter(peak==el & term=='I' & estimate >0) %>% distinct(aa))$aa
        good_aa <- intersect(good_aa, el_aa[[el]])
        i0s <- (dd %>% filter(aa %in% good_aa & term=='I' & peak=='0'))$estimate
        i1s <- (dd %>% filter(aa %in% good_aa & term=='I' & peak==el) %>% mutate(x=estimate/n))$x
        i0s <- i0s[i1s>0]
        i1s <- i1s[i1s>0]
        rs <- i1s/i0s
        r <- sum(i1s[abs(rs-median(rs))<1.5*mad(rs)])/sum(i0s[abs(rs-median(rs))<1.5*mad(rs)])
        res <- res %>% bind_rows(data.frame(element=el, gamma=r, logI = log2(sum(i0s))))
      }
      res$scan = dd$scan[[1]]
      res$rt = dd$rt[[1]]
      res
    }) %>% ungroup() -> any_aa
  any_aa %>%
    filter(rt>900, rt<3800) %>%
    filter(gamma<0.02 & gamma>1e-5) %>%
    ggplot(aes(x=logI, y=gamma)) + geom_point()+geom_smooth() + theme_bw() + facet_wrap(~element, scales = 'free')


  result %>%
    select(term, estimate, peak, aa, rt, scan, n) %>%
    group_by(scan, rt, aa) %>%
    do({
      dd <- .
      res <- data.frame()
      elements <- setdiff(unique(dd$peak), '0')
      for(el in elements){
        if(nrow(dd %>% filter(term=='I' & peak==el))>0){
          i0 <- sum((dd %>% filter(term=='I' & peak=='0'))$estimate)
          i1 <-sum((dd %>% filter(term=='I' & peak==el) %>% mutate(x=estimate/n))$x)
          res <- res %>% bind_rows(data.frame(element=el, gamma=i1/i0, logI = log2(i0)))
        }
      }
      res$scan = dd$scan[[1]]
      res$rt = dd$rt[[1]]
      res$aa = dd$aa[[1]]
      res
    }) %>% ungroup() -> all_aa
  all_aa %>%
    bind_rows(any_aa %>% mutate(aa='any')) %>%
    filter(aa %in% c('L','P','V','any')) %>%
    filter(rt>900, rt<3800) %>%
    filter(gamma<0.02 & gamma>1e-5) %>%
    ggplot(aes(x=rt, y=gamma)) + geom_point()+geom_smooth() + theme_bw() + facet_grid(element~aa, scales = 'free')

  result %>%
    select(term, estimate, peak, aa, rt, scan, n) %>%
    group_by(scan, rt, aa) %>%
    do({
      dd <- .
      res <- data.frame()
      elements <- setdiff(unique(dd$peak), '0')
      if('C' %in% elements){
        for(el in elements){
          i1 <- (dd %>% filter(term=='I' & peak==el) %>% mutate(estimate=estimate/n))$estimate
          i0 <- (dd %>% filter(term=='I' & peak=='0'))$estimate
          iC <- (dd %>% filter(term=='I' & peak=='C') %>% mutate(estimate=estimate/n))$estimate
          res <- res %>% bind_rows(data.frame(element=el, gamma=i1/i0, gammaC = i1/iC, logI = log2(i0)))
        }
        res$scan = dd$scan[[1]]
        res$rt = dd$rt[[1]]
        res$aa = dd$aa[[1]]
      }
      res
    }) %>% ungroup() -> all_aa
  all_aa %>%
    filter(aa %in% c('L','P','V','H')) %>%
    filter(rt>900, rt<3800) %>%
    filter(gamma<0.02 & gamma>1e-5) %>%
    ggplot(aes(x=rt, y=gammaC*0.0110)) + geom_point()+geom_smooth() + theme_bw() + facet_grid(element~aa, scales = 'free')
  all_aa %>%
    filter(aa %in% c('L','P','V','H')) %>%
    filter(rt>900, rt<3800) %>%
    filter(gamma<0.02 & gamma>1e-5) %>%
    ggplot(aes(x=rt, y=gamma)) + geom_point()+geom_smooth() + theme_bw() + facet_grid(element~aa, scales = 'free')

  any_aa %>%
    filter(rt>900, rt<3800) %>%
    filter(gamma<0.02 & gamma>1e-5) %>%
    group_by(element) %>%
    summarise(mean=mean(gamma), median=median(gamma), sd=sd(gamma), se=sd(gamma)/sqrt(n()-1), meanI=mean(logI), medianI=median(logI)) %>%
    mutate(file=file)
}

fpath <- "Fusion/Methods/20170617/"
files <- data_frame(path=fpath, file=list.files(fpath, pattern="^2017.*.mzML"))

files %>%
  group_by(file) %>%
  do({
    dd <- analyze_immoniums(file=file.path(.$path[[1]],.$file[[1]]))
    write_csv(dd, sub(".mzML","_gaussian.csv",file.path(.$path[[1]],.$file[[1]]),))
    dd
  }) -> all_results

all_results %>%
  mutate(sample=sub(".*tMS2_(\\d+uscans)_(\\d+us)_([^_]+)_(\\d+K)_(.*).mzML","\\4_\\5", file)) %>%
  mutate(sample=sub(".*EColi_","", sample)) %>%
  mutate(sample=sub("o1","O1", sample)) %>%
  # mutate(mass=sub(".*_(\\d+)ng.mzML","\\1", file)) %>%
  # mutate(mass=as.numeric(mass)) %>%
  #filter(!grepl("mzML", sample)) %>%
  # filter(!grepl("20170419_F_13_HeLa_Imm", file)) %>%
  ggplot(aes(x=(sample), y=median, ymin = median-se, ymax=median+se, color=meanI))+geom_point() +geom_errorbar(size=0.1)+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))+ facet_grid(element~., scales='free')
all_results %>%
  mutate(sample=sub(".*tMS2_(\\d+us)_(\\d+K)_(.*).mzML","\\1_\\2_\\3", file)) %>%
  mutate(sample=sub(".*EColi_","", sample)) %>%
  mutate(sample=sub("o1","O1", sample)) %>%
  # mutate(mass=sub(".*_(\\d+)ng.mzML","\\1", file)) %>%
  # mutate(mass=as.numeric(mass)) %>%
  #filter(!grepl("mzML", sample)) %>%
  # filter(!grepl("20170419_F_13_HeLa_Imm", file)) %>%
  ggplot(aes(x=(meanI), y=median, ymin = median-se, ymax=median+se, color=sample)) + geom_point() +geom_errorbar()+theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))+ facet_grid(element~., scales='free')

write_csv(result, "/Volumes/Seagate 2T/Immonium/Fusion/HeLa/Conc/20170408_F_04_HeLa_Imm_200ng.peakfitting.csv")
result %>%
  filter(peak=="Mono" & term=='I') %>%
  select(i0=estimate, ion, scan) -> i0
#i0$C <- apply(i0, 1, function(x)immoniumIons[[x["ion"]]]['C'])
#i0$N <- apply(i0, 1, function(x)immoniumIons[[x["ion"]]]['N'])
#i0$H <- apply(i0, 1, function(x)immoniumIons[[x["ion"]]]['H'])
#i0$O <- apply(i0, 1, function(x)immoniumIons[[x["ion"]]]['O'])
#i0$O[is.na(i0$O)] <- 0

result %>%
  filter(term=='I') %>%
  select(i=estimate, ion, peak, scan, rt) %>%
  left_join(i0, by=c('ion', 'scan'))-> data

data$n <- apply(data, 1, function(x)immoniumIons[[x["ion"]]][x["peak"]])

data %>%
  filter(!is.na(n)) %>%
  mutate(gamma=i/i0/n) %>%
  mutate(logI = log2(i0)) %>%
  filter(gamma<0.015 & rt>900 & rt < 3800) %>%
  ggplot(aes(x=logI, y=gamma)) + geom_point() + geom_smooth() + facet_grid(peak~ion, scales = 'free_y')
