###########################################
#'@export
###########################################
get.GLD <- function(path='.', loci=NULL, exclude = 'ALL', ...){
  #Some initial parameter setup and autodetection for the populations and loci
  pops <- setdiff(dir(path) , dir(path, pattern='\\.'))
  pops <- setdiff(pops, exclude)
  name <- unlist(strsplit(dir(file.path(path, pops[1]), pattern='.ga2l$')[1], paste0('_',pops[1],'_|\\.ga2l')))[1]
  if(is.null(loci)){
    for (pop in pops){
      freq.files <- dir(file.path(path, pop), pattern='.ga2l$')
      for (freq.file in freq.files){
        loci <- c(loci,unlist(strsplit(freq.file, paste0('_',pop,'_|\\.ga2l')))[2])
        #lines <- readLines(file(file.path(path, pop, freq.file), open='r'))
      }
    }
    loci <- sort(unique(loci))
  }
  
  #
  popFreq <- list()
  for (locus in loci){
    LRT <- c()
    PRS <- c()
    for (pop in pops){
      filename <- paste(name, pop, paste0(locus, '.ga2l'),sep='_')
      f <- file(file.path(path, pop, filename), open='r')
      temp <- readLines(f)
      close(f)
      LRT <- c(LRT, as.numeric(unlist(strsplit(temp[1], 'lnL ratio test statistic | p-value on |\\. df: '))[4]))
      PRS <- c(PRS, as.numeric(unlist(strsplit(temp[2], '  \\(observed p-value is quantile | of | simulations\\)'))[2]))
    }
    if(is.null(popFreq$LRT)){
      popFreq$LRT <- data.frame(LRT, row.names=pops)
      popFreq$PRS <- data.frame(PRS, row.names=pops)
    } else {
      popFreq$LRT <- cbind(popFreq$LRT, LRT)
      popFreq$PRS <- cbind(popFreq$PRS, PRS)
    }
  }
  closeAllConnections()
  colnames(popFreq$LRT) <- loci
  colnames(popFreq$PRS) <- loci
  return(popFreq)
}
