###########################################
#'@export
###########################################
get.outliers <- function(tmp, exclude=NULL){
  far <- compare.further(tmp, exclude=exclude)
  for (n in names(far)){
    far[[n]] <- far[[n]][apply(far[[n]], 1, any),]
  }
  return(far)
}

###########################################
#'@export
###########################################
unlist.outliers <- function(far){
  outl = c()
  for (n in names(far)){
    outl <- c(outl, unlist(lapply(apply(far[[n]], 1, which), paste, collapse=' ')))
  }
  return(cbind(names(outl), outl))
}

###########################################
#'@export
###########################################
count.outliers <- function(far){
  df <- data.frame()
  for (n in names(far)){
    temp <- apply(far[[n]], 2, sum)
    partDF <- data.frame(names(temp),temp)
    names(partDF) <- c('region', n)
    if ('region' %in% names(df)){
      df <- merge(df, partDF,all=TRUE)
    } else {
      df <- partDF
    }
  }
  return(df)
}
