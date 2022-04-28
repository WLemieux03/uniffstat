###########################################
#'@export
###########################################
EW.calc <- function(path='.', method=c('Slatkin', 'Ewens'), side=c('low', 'high', 'both'), loci=NULL, exclude = 'ALL', ...){
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
  method <- match.arg(method)
  side <- match.arg(side)
  #
  EW <- NULL
  for (pop in pops){
    poptemp <- c()
    for (locus in loci){
      filename <- paste(name, pop, paste0(locus, '.ewres'),sep='_')
      f <- file(file.path(path, pop, filename), open='r')
      temp <- readLines(f)
      close(f)
      neut <- as.numeric(unlist(strsplit(temp[grep(method, temp)+1], ' '))[c(7,9,11,13)])
      
      if (side=='low'){
        res <- (neut[1]+c(0,neut[3])[(neut[1]!=0)+1])
      } else if (side=='high'){
        res <- (neut[2]+c(0,neut[4])[(neut[2]!=0)+1])
      } else if (side=='both'){
        res <- paste((neut[1]+c(0,neut[3])[(neut[1]!=0)+1]), (neut[2]+c(0,neut[4])[(neut[2]!=0)+1]), sep = '/')
      }
      poptemp <- c(poptemp, res)
    }
    EW <- cbind(EW, poptemp)
  }
  colnames(EW) <- pops
  rownames(EW) <- loci
  closeAllConnections()
  return(EW)
}
