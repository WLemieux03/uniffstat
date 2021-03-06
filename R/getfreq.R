###########################################
#'@export
###########################################
get.freq <- function(path='.', loci=NULL, exclude = 'ALL', ...){
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
  }
  
  #
  popFreq <- list()
  for (locus in loci){
    locusDF <- data.frame()
    for (pop in pops){
      filename <- paste(name, pop, paste0(locus, '.freq'),sep='_')
      f <- file(file.path(path, pop, filename), open='r')
      temp <- readLines(f)
      close(f)
      alleles <- c()
      freq <- c()
      for (t in temp[grep(paste0(locus, '\\*'), temp)]){
        all_freq <- unlist(strsplit(t, '\t'))
        alleles <- c(alleles, all_freq[1])
        freq <- c(freq, as.numeric(all_freq[2]))
      }
      partDF <- data.frame(alleles,freq)
      names(partDF) <- c('alleles', pop)
      if('alleles' %in% names(locusDF)){
        locusDF <- merge(locusDF, partDF, all=T)
      } else {
        locusDF <- partDF
      }
    }
    row.names(locusDF) <- locusDF$alleles
    locusDF[is.na(locusDF)]=0
    popFreq[[locus]] <- as.matrix(locusDF[,-1])
  }
  closeAllConnections()
  return(popFreq)
}
