if(F){
  library(isoms)

  fpath <- "~/tmp/Immonium"
  files <- data_frame(path=fpath, file=list.files(fpath, pattern="^2017.*fullP.*.mzML"))
  ff = apply(files, 1, function(x){file.path(x['path'],x['file'])})
  ncpu=8
  result <- bind_rows(mclapply(
      ff, function(ii){
        outfile <- sub(".mzML","_gaussian.csv", ii)
        if(!file.exists(outfile)){
          xxx <- analyze_immoniums(ii, ions=immoniumIons, width=0.0015, fixSigma=T)
          write_csv(xxx, outfile)
        }else
          xxx <- read_csv(outfile)
        xxx %>% mutate(file=outfile)
      }, mc.cores=ncpu, mc.preschedule=TRUE))


  analyze_immoniums(file="~/tmp/Immonium/20170626_05_tMS2only_10uscans_20us_5e4_60K_10SF_200-1000_EColi_400.mzML",
                    ions = immoniumIons) -> xxx

  xxx %>%
    filter(R2>0.95) %>%
    ggplot(aes(x=isoratio, fill=peak)) + geom_density(alpha=0.3) + facet_grid(peak~ion, scales = 'free')
  xxx %>%
    filter(sigma>0&R2>0.95) %>%
    ggplot(aes(x=mu, y=sigma, color=ion, shape=peak)) + geom_point()# + facet_grid(peak~ion, scales = 'free')

  xxx %>%
    filter(peak %in%c('13C','15N','2H','18O')) %>%
    filter(ion %in% c('L','P','V','T','Pyr','F','H')) %>%
    filter(isoratio >0 & abs(masserror)<1e-4) %>%
    ggplot(aes(x=(isoratio/n), fill=ion)) + geom_density(alpha=0.3) + facet_grid(.~peak, scales = 'free')
  xxx %>%
    group_by(seqNum) %>%
    summarize(n=length(unique(ion))) %>% ggplot(aes(x=seqNum, y=n)) + geom_point()
  xxx %>%
    group_by(seqNum, rt) %>%
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
        good_aa <- (dd %>% filter(peak==el & I >0) %>% distinct(ion))$ion
#        good_aa <- intersect(good_aa, el_aa[[el]])
        i0s <- (dd %>% filter(ion %in% good_aa & peak=='0'))$I
        i1s <- (dd %>% filter(ion %in% good_aa & peak==el) %>% mutate(x=I/n))$x
        rs <- i1s/i0s
        i0s <- i0s[rs<1]
        i1s <- i1s[rs<1]
        rs <- rs[rs<1]
        rs_good <- (abs(rs-median(rs))<1.5*mad(rs))
        r <- sum(i1s[rs_good])/sum(i0s[rs_good])
        if(!is.na(r))
          res <- res %>% bind_rows(data.frame(peak=el, gamma=r, logI = log2(sum(i0s))))
      }
      if(nrow(res)>0){
        res$seqNum = dd$seqNum[[1]]
        res$rt = dd$rt[[1]]
        res$ion = 'any'
      }
      res
    }) %>% ungroup() -> any_aa
  any_aa %>%
    filter(peak %in%c('15N')) %>%
    filter(gamma<0.01) %>%
#    filter(ion %in% c('L','P','V','T','Pyr','F','H')) %>%
#    filter(isoratio >0 & abs(masserror)<1e-4) %>%
    ggplot(aes(x=((gamma)), fill=peak)) + geom_density(bins = 200)
  xxx %>%
    filter(peak %in%c('15N')) %>%
    mutate(gamma=isoratio/n) %>%
    filter(gamma<0.01) %>%
    #    filter(ion %in% c('L','P','V','T','Pyr','F','H')) %>%
    #    filter(isoratio >0 & abs(masserror)<1e-4) %>%
    ggplot(aes(x=gamma, fill=ion)) + geom_density(bins = 200)

}
