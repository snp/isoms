analyze_immoniums <- function(file, width=0.001, ions=immoniumIons, fixSigma=T){
  message(sprintf("Reading file [%s]", file))
  msrun <- openMSfile(file, backend = "Ramp")
  hd <- header(msrun)
  immscans <- which(hd$collisionEnergy>45 & (hd$msLevel>1))
  message(sprintf("\t%d immonium scans found", length(immscans)))
  result <- data.frame()
  message("Processing:")

  hd %>%
    filter(seqNum %in% immscans) %>%
    group_by(seqNum) %>%
    do({
      ires <- data.frame()
      dd = .
      ii <- dd$seqNum[[1]]
      ss <- peaks(msrun, ii)
      ss_range <- diff(range(ss[,1]))
      for(ion_ in names(ions)){
        mz = monoMass(ions[[ion_]])
        peaks <- c("13C")
        if("N" %in% names(ions[[ion_]]) & ions[[ion_]]["N"] > 0)
          peaks <- c(peaks, "15N")
        if("H" %in% names(ions[[ion_]]) & ions[[ion_]]["H"] > 0)
          peaks <- c(peaks, "2H")
        if("O" %in% names(ions[[ion_]]) & ions[[ion_]]["O"] > 0)
          peaks <- c(peaks, "18O")
        if(ss_range>5){
          cres <- get_isopeaks(ss, mz, width=width, npoint=6, fixSigma=fixSigma, peaks=peaks)
        }else{
          cres <- get_isopeaks_nomono(ss, mz, width=width, npoint=6, fixSigma=fixSigma)
        }
        if(nrow(cres)>1)
          ires <- ires %>%
            bind_rows(
              cres %>%
                mutate(
                  rt = dd$retentionTime[[1]],
                  tic = dd$totIonCurrent[[1]],
                  bpI = dd$basePeakIntensity[[1]],
                  bpMZ = dd$basePeakMZ[[1]],
                  ion=ion_)
            )
      }
      ires$n <- apply(ires, 1, function(x)ions[[x["ion"]]][sub("\\d+","",x["peak"])])
      ires
    }) -> result
  result
}

analyze_immoniums2 <- function(file, width=0.0015, ions=immoniumIons){
  message(sprintf("Reading file [%s]", file))
  msrun <- openMSfile(file, backend = "Ramp")
  hd <- header(msrun)
  immscans <- which(hd$collisionEnergy>49 & (hd$msLevel>1))
  message(sprintf("\t%d immonium scans found", length(immscans)))
  result <- data.frame()
  message("Processing:")

  hd %>%
    filter(seqNum %in% immscans) %>%
    group_by(seqNum) %>%
    do({
      ires <- data.frame()
      dd = .
      ii <- dd$seqNum[[1]]
      ss <- peaks(msrun, ii)

      ss_range <- diff(range(ss[,1]))
      for(ion_ in names(ions)){
        mz = monoMass(ions[[ion_]])
        if(ss_range>5){
          cres <- get_isopeaks2(ss, mz, width=width, npoint=6)
        }
        if(nrow(cres)>0)
          ires <- ires %>%
            bind_rows(
              cres %>%
                mutate(
                  rt = dd$retentionTime[[1]],
                  tic = dd$totIonCurrent[[1]],
                  bpI = dd$basePeakIntensity[[1]],
                  bpMZ = dd$basePeakMZ[[1]],
                  ion=ion_)
            )
      }
      ires$element <- toupper(sub("r","",ires$term))
      ires$n <- apply(ires, 1, function(x)ions[[x["ion"]]][x["element"]])
      ires
    }) -> result
  result
}

analyze_immoniums3 <- function(file, width=0.0015, ions=immoniumIons){
  message(sprintf("Reading file [%s]", file))
  msrun <- openMSfile(file)
  hd <- header(msrun)
  immscans <- which(hd$collisionEnergy>49 & (hd$msLevel>1))
  message(sprintf("\t%d immonium scans found", length(immscans)))
  result <- data.frame()
  message("Processing:")
  result <- bind_rows(mclapply(immscans, function(ii){
    msrun <- openMSfile(file)
    ss <- peaks(msrun, ii)
    ss_range <- diff(range(ss[,1]))
    for(ion_ in names(ions)){
      mz = monoMass(ions[[ion_]])
      if(ss_range>5){
        cres <- get_isopeaks2(ss, mz, width=width, npoint=6)
      }
      if(nrow(cres)>0)
        ires <- ires %>%
          bind_rows(
            cres %>%
              mutate(ion=ion_)
          )
    }
    ires$element <- toupper(sub("r","",ires$term))
    ires$n <- apply(ires, 1, function(x)ions[[x["ion"]]][x["element"]])
    ires
  }, mc.cores=detectCores(), mc.preschedule=TRUE))

  hd %>%
    filter(seqNum %in% immscans) %>%
    group_by(seqNum) %>%
    do({
      ires <- data.frame()
      dd = .
      ii <- dd$seqNum[[1]]
      message(ii)
      ss <- peaks(msrun, ii)

      ss_range <- diff(range(ss[,1]))
      for(ion_ in names(ions)){
        mz = monoMass(ions[[ion_]])
        if(ss_range>5){
          cres <- get_isopeaks2(ss, mz, width=width, npoint=6)
        }
        if(nrow(cres)>0)
          ires <- ires %>%
            bind_rows(
              cres %>%
                mutate(
                  rt = dd$retentionTime[[1]],
                  #                  scan = ii,
                  ion=ion_)
            )
      }
      ires$element <- toupper(sub("r","",ires$term))
      ires$n <- apply(ires, 1, function(x)ions[[x["ion"]]][x["element"]])
      ires
    }) -> result
  result
}

