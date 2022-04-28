get.fullHaplo <- function(path='.',locus='A~B~C~DRB1~DQB1~DPB1',exclude = 'ALL', ...){
  #Some initial parameter setup and autodetection for the populations and loci
  pops <- setdiff(dir(path) , dir(path, pattern='\\.'))
  pops <- setdiff(pops, exclude)
  name <- unlist(strsplit(dir(file.path(path, pops[1]), pattern='.freq$')[1], paste0('_',pops[1],'_|\\.freq')))[1]
  #
  locusDF <- data.frame()
  for (pop in pops){
    filename <- paste(name, pop, paste0(locus, '.freq'),sep='_')
    partDF <- read.table(file.path(path, pop, filename), skip=6)
    names(partDF) <- c('alleles', pop)
    if('alleles' %in% names(locusDF)){
      locusDF <- merge(locusDF, partDF, all=T)
    } else {
      locusDF <- partDF
    }
    
    
  }
  row.names(locusDF) <- locusDF$alleles
  locusDF[is.na(locusDF)]=0
  locusDF <- locusDF[apply(locusDF[,-1], 1, function(x){any(x!=0)}),]
  attributes(locusDF)$sum <- locusDF['sum:',-1]
  attributes(locusDF)$class <- c(attributes(locusDF)$class, 'locusDF')
  locusDF <- locusDF[setdiff(rownames(locusDF), 'sum:'),]
  return(locusDF)
}
