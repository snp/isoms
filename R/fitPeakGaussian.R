
fitPeak <- function(x, y, startPars=c('mu' = 70.065, 'I' = 5e3, 'sigma' = 3e-4), maxAttempts=1000, sigma=0){
  fctStr <- function(x, mu, I, sigma) { I/sigma*exp(-1/2*(x-mu)^2/sigma^2)} # "I/sigma*exp(-1/2*(x-mu)^2/sigma^2)"
  fitFct <- as.formula("y ~ fctStr(x, mu, I, sigma)")
  if(sigma>0){
    fctStr <- function(x, mu, I) {I/sigma*exp(-1/2*(x-mu)^2/sigma^2)}
    startPars <- startPars[c('mu','I')]
    fitFct <- as.formula("y ~ fctStr(x, mu, I)")
  }
#  fitFct <- as.formula(paste("y ~", fctStr))
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
