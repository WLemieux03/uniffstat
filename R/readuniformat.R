
get.AAM <- function(x, f){
  require(reshape2)
  if(x=='@'){
    domain <- names(f[f!=0])
    haplo <- expand.grid(domain, domain)
    colnames(haplo) <- c('A1','A2')
  }else{
    domain <- sort(unique(unlist(strsplit(x, '\\||,'))))
    haplo <- unlist(strsplit(x, '\\|'))
    haplo <- colsplit(haplo, ',',c('A1','A2'))
  }
  haplo <- cbind(haplo, sapply(domain,function(y){apply(haplo, 1, function(x){y %in% x})}, simplify=F))
  tmp <- t(t(haplo[,domain])*f[domain])
  tmp[is.na(tmp)] <- 0
  tmp[!haplo[domain]] <- 1
  if (length(domain)==1){
    mult <- c(2)
  } else {
    mult <- (3-apply(haplo[,domain],1,sum))
  }
  haplo <- cbind(haplo, freq=apply(tmp,1,prod)^mult, mult=mult)
  res <- as.data.frame(t(sapply(domain, function(y){sum(haplo[haplo[,y],'freq']*haplo[haplo[,y],'mult'])/sum(haplo[,'freq'])})))
  return(res)
}


UNFMTtoTBL <- function(temp, f){
  require(foreach)
  x <- foreach(i=1:length(temp)) %do% {tmp = get.AAM(temp[i], f); rownames(tmp)=i;tmp}
  res <- NULL
  for (i in 1:length(x)){
    if(is.null(res)){
      res <- cbind(rn=row.names(x[[i]]), x[[i]])
    }else{
      res <- merge(res, cbind(rn=row.names(x[[i]]), x[[i]]), all=TRUE)
    }
  }
  res[is.na(res)] <- 0
  row.names(res) <- res[,1]
  res <- res[,-1]
  res <- res[,colSums(res)!=0]
  return(res[order(as.integer(rownames(res))),sort(colnames(res))])
}

###########################################
#'@export
###########################################
read.uniformat <- function(path='.', freq=NULL, exclude='ALL', coll=F, ...){
  fs <- dir(path, pattern='\\.unif')
  IDs <- NULL
  rn <- NULL
  data <- list()
  for (f in fs){
    tmp <- unlist(strsplit(f, '_|\\.'))
    pop <- tmp[3]
    if(pop %in% exclude){next}
    loci <- unlist(strsplit(tmp[4], '~'))
    dataf <- read.table(file.path(path, f))
    IDs <- c(IDs, dataf[,1])
    rn <- c(rn, rep(pop, length(dataf[,1])))
    for (l in 1:length(loci)){
      temp <- UNFMTtoTBL(dataf[,l+1], freq[[loci[l]]][,pop])
      if (is.null(data[[loci[l]]])){
        data[[loci[l]]] <- cbind(id=dataf[,1],temp)
      } else {
        data[[loci[l]]] <- merge(data[[loci[l]]], cbind(id=dataf[,1],temp), all=T)
      }
      
    }
  }
  data <- lapply(data, function(x){x[,-1]})
  data <- lapply(data, function(x){x[is.na(x)] <- 0;x})
  if(!coll){
    return(list(id=IDs, rn=rn, data=data))
  } else {
    return(collate(list(id=IDs, rn=rn, data=data)))
  }
  
}


###########################################
#'@export
###########################################
collate <- function(data){
  dos <- data$data[[names(data$data)[1]]]
  for (n in names(data$data)[-1]){
    dos <- cbind(dos,data$data[[n]])
  }
  dos
}