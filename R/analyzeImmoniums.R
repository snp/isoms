analyze_immoniums <- function(file, width=0.001, ions=immoniumIons, fixSigma=T){
  message(sprintf("Reading file [%s]", file))
  msrun <- openMSfile(file)
  hd <- header(msrun)
  immscans <- which(hd$collisionEnergy>49 & (hd$msLevel>1))
  message(sprintf("\t%d immonium scans found", length(immscans)))
  result <- data.frame()
  message("Processing:")
#   require(multidplyr)
#   cl <- create_cluster()
#   cluster_library(cl, "isoms")
#   cluster_copy(cl, file)
#   cluster_copy(cl, isotopes)
#   cluster_copy(cl, width)
#   cluster_copy(cl, ions)
#   cluster_copy(cl, immoniumIons)
#   cluster_assign_expr(cl, "msrun", quote(openMSfile(file)))
# #  cluster_copy(cl, immscans)
#   set_default_cluster(cl)

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
          cres <- get_isopeaks(ss, mz, width=width, npoint=6, fixSigma=fixSigma)
        }else{
          cres <- get_isopeaks_nomono(ss, mz, width=width, npoint=6, fixSigma=fixSigma)
        }
        if(nrow(cres)>1)
          ires <- ires %>%
            bind_rows(
              cres %>%
                mutate(
                  rt = dd$retentionTime[[1]],
#                  scan = ii,
                  ion=ion_)
            )
      }
      ires$n <- apply(ires, 1, function(x)ions[[x["ion"]]][sub("\\d+","",x["peak"])])
      ires
    }) -> result
  result
}
