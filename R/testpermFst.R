###########################################
#'@export
###########################################
test.perm <- function(data, pop, nperm=1000, ...){
  npop=length(unique(pop))
  perm.stat <- matrix(nrow=npop+1, ncol=nperm)
  pw <- matrix(1,nrow=npop, ncol=npop)
  
  fs <- fs.dosage(data, pop, ...)
  perm.stat[,nperm] <- fs$Fs['Fst',]
  rownames(perm.stat) <- names(fs$Fs['Fst',])
  rownames(pw) <- colnames(pw) <- names(fs$Fs['Fst',])[1:npop]
  for (i in 1:(nperm - 1)) {
    fsp <- fs.dosage(data, sample(pop), ...)
    perm.stat[,i] <- fsp$Fs['Fst',]
    pw <- pw+(fsp$Fst2x2>=fs$Fst2x2)
  }
  res <- apply(perm.stat, 1, function(x){sum(x >= x[nperm])/nperm})
  names(res) <- rownames(perm.stat)
  pw <- pw/nperm
  list(p.val=res, pairwise=pw, Fs=fs)
}

###########################################
#'@export
###########################################
matching <- function (dos) {
  if (!is.matrix(dos)) {
    if (class(dos)[[1]] == "bed.matrix") 
      dos <- gaston::as.matrix(dos)
    else dos <- as.matrix(dos)
  }
  lims <- range(dos, na.rm = TRUE)
  if ((lims[2] > 2) | (lims[1] < 0)) 
    stop("input dosage matrix should contains only 0, 1 and 2s")
  if (sum(is.na(dos)) > 0) {
    na <- matrix(rep(1, prod(dim(dos))), ncol = ncol(dos))
    ina <- which(is.na(dos))
    na[ina] <- 0
    dos[ina] <- 1
    Mij <- 1/2 * (1 + 1/tcrossprod(na) * tcrossprod(dos - 
                                                      1))
  }
  else {
    nl <- dim(dos)[2]
    Mij <- 1/2 * (1 + tcrossprod(dos - 1)/nl)
  }
  Mij
}
