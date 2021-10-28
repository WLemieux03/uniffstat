###########################################
#'@export
###########################################
merge.freq <- function(X, Y){
  loci <- unique(c(names(X),names(Y)))
  popFreq <- list()
  for (locus in loci){
    x <- cbind(data.frame(alleles=rownames(X[[locus]])), X[[locus]])
    y <- cbind(data.frame(alleles=rownames(Y[[locus]])), Y[[locus]]) 
     
    locusDF <- merge(x, y, all=T)
    row.names(locusDF) <- locusDF$alleles
    locusDF[is.na(locusDF)]=0
    popFreq[[locus]] <- as.matrix(locusDF[,-1])
  }
  return(popFreq)
}
