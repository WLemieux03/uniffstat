###########################################
#'@export
###########################################
rfreq <- function(alleleList=sampleList, loci=names(alleleList), pops=c('01','02')){
  popFreq <- list()
  for (locus in loci){
    locusDF <- data.frame()
    for (pop in pops){
      alleles <- sort(sample(alleleList[[locus]], sample(1:length(alleleList[[locus]]), 1)))
      freq <- sapply(rnorm(length(alleles), 0.3,0.3), function(x){max(x,0)})
      freq <- freq/sum(freq)
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
  return(clean.freq(popFreq))
}

