###########################################
#'@export
###########################################
loci.dist<-function(path='.',diploid=TRUE,method="Dch", loci=c('A', 'B', 'C', 'DPB1', 'DQB1', 'DRB1'), exclude='ALL', ...){
  match.call()
  by.loci <- list()
  for (locus in loci){
    by.loci[[locus]] <- genet.dist(path=path, diploid=diploid, method=method, loci=locus, exclude=exclude)
  }
  global <- genet.dist(path=path, diploid=diploid, method=method, loci=loci, exclude=exclude)
  results <- list(dist.by.loci=by.loci, dist.global=global)
  return(results)
}