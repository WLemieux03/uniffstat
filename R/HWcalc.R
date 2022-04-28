###########################################
#'@export
###########################################
HW.calc <- function(path='.', loci=NULL, exclude = 'ALL', ...){
  #Some initial parameter setup and autodetection for the populations and loci
  pops <- setdiff(dir(path) , dir(path, pattern='\\.'))
  pops <- setdiff(pops, exclude)
  name <- unlist(strsplit(dir(file.path(path, pops[1]), pattern='.freq$')[1], paste0('_',pops[1],'_|\\.freq')))[1]
  if(is.null(loci)){
    for (pop in pops){
      freq.files <- dir(file.path(path, pop), pattern='.freq$')
      for (freq.file in freq.files){
        loci <- c(loci,unlist(strsplit(freq.file, paste0('_',pop,'_|\\.freq')))[2])
        #lines <- readLines(file(file.path(path, pop, freq.file), open='r'))
      }
    }
    loci <- sort(unique(loci))
    loci <- loci[-grep('~', loci)]
  }
  
  #
  HW <- NULL
  for (pop in pops){
    poptemp <- c()
    for (locus in loci){
      filename <- paste(name, pop, paste0(locus, '.freq'),sep='_')
      f <- file(file.path(path, pop, filename), open='r')
      temp <- readLines(f)
      close(f)
      cgc <- as.numeric(unlist(strsplit(temp[grep('max lnL', temp)], '\t'))[2])
      
      filename <- paste(name, pop, paste0(locus, '.freqi'),sep='_')
      f <- file(file.path(path, pop, filename), open='r')
      temp <- readLines(f)
      close(f)
      cgci <- as.numeric(unlist(strsplit(temp[grep('max lnL', temp)], '\t'))[2])
      
     poptemp <- c(poptemp, abs(cgc-cgci))
    }
    HW <- cbind(HW, poptemp)
  }
  colnames(HW) <- pops
  rownames(HW) <- loci
  closeAllConnections()
  return(HW)
}
