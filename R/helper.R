further <- function(x, type=c('median', 'quartile')){
  type <- match.arg(type)
  iqr <- IQR(x)
  med <- median(x)
  q <- quantile(x, c(.25, .75))
  if (type=='median'){
    return((x < (med-1.5*iqr)) | (x > (med+1.5*iqr)))
  } else if (type=='quartile'){
    return((x < (q[1]-1.5*iqr)) | (x > (q[2]+1.5*iqr)))
  }
}

compare.further <- function(tmp, exclude=NULL){
  res <- list()
  for (n in names(tmp)){
    regions <- setdiff(colnames(tmp[[n]]), exclude)
    res[[n]] <- t(apply(tmp[[n]][,regions],1,further))
  }
  return(res)
}

remove.rare <- function(dfl, threshold=0.01, disr=NULL){
  for (n in names(dfl)){
    mask <- setdiff(colnames(dfl[[n]]), disr)
    dfl[[n]] <- dfl[[n]][apply(dfl[[n]][,mask], 1, function(x){any(x>=threshold)}),]
  }
  return(dfl)
}

which.diff <- function(dfl, x, exclude=NULL){
  DF <- data.frame()
  FF <- data.frame()
  far <- get.outliers(dfl, exclude=exclude)
  for (n in names(dfl)){
    if(length(x)>1){
      mask <- names(apply(far[[n]][,x], 1, any))
    } else if(length(x)==1) {
      mask <- names(far[[n]][,x])
    } else {
      mask <- rownames(far[[n]])
    }
    tmp <- dfl[[n]][mask,]
    ftmp <- far[[n]][mask,]
    if (dim(DF)[1]==0){
      DF <- tmp
      FF <- ftmp
    } else {
      DF <- rbind(DF, tmp)
      FF <- rbind(FF, ftmp)
    }
  }
  output <- list(freq=DF, status=FF)
  return(output)
}

list.merge <- function(X, Y, exclude=NULL){
  popFreq <- list()
  for (locus in union(names(X), names(Y))){
    
    locusDF <- cbind(data.frame(allele=rownames(X[[locus]])),X[[locus]][, setdiff(colnames(X[[locus]]), exclude)])
    if (dim(locusDF)[2]<3){colnames(locusDF) <- c("allele", setdiff(colnames(X[[locus]]), exclude))}
    partDF <- cbind(data.frame(allele=rownames(Y[[locus]])),Y[[locus]][, setdiff(colnames(Y[[locus]]), exclude)])
    if (dim(partDF)[2]<3){colnames(partDF) <- c("allele", setdiff(colnames(Y[[locus]]), exclude))}
    locusDF <- merge(locusDF, partDF, all=T)
    
    row.names(locusDF) <- locusDF$allele
    locusDF[is.na(locusDF)]=0
    popFreq[[locus]] <- as.matrix(locusDF[,-1])
  }
  return(popFreq)
}
