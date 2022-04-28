pasteUnif <- function(ind, loci=c('A', 'B', 'C', 'DRB1', 'DQB1', 'DPB1'), dict=NULL){
  require(tidyr)
  require(curl)
  ustring <- ind[["IND"]]
  if (is.null(dict)){
    dict <- list()
  }
  for (locus in loci){
    if (any(is.na(ind[paste0(locus, c("_1","_2"))]))){
      ustring <- c(ustring, "@")
      next
    }
    A1 <- paste(locus, ind[paste0(locus, "_1")], sep="*")
    A2 <- paste(locus, ind[paste0(locus, "_2")], sep="*")
    if (grepl("\\d{2,4}:[[:alpha:]]{2,5}|UNK\\d{4}", A1)){
      if (A1 %in% names(dict)){
        A1 <- dict[[A1]]
      }else{
        con <- curl(paste0("https://hml.nmdp.org/mac/api/decode?typing=", A1))
        temp <- tryCatch({unlist(strsplit(readLines(con), "/"))}, 
                         error=function(cond){message(paste("MAC", A1, "does not exist")); message(cond); return(NULL)})
        close(con)
        dict[[sear]] <- temp
        A1 <- temp
      }
    }
    if (grepl("\\d{2,4}:[[:alpha:]]{2,5}|UNK\\d{4}", A2)){
      if (A2 %in% names(dict)){
        A2 <- dict[[A2]]
      }else{
        con <- curl(paste0("https://hml.nmdp.org/mac/api/decode?typing=", A2))
        temp <- tryCatch({unlist(strsplit(readLines(con), "/"))}, 
                         error=function(cond){message(paste("MAC", A2, "does not exist")); message(cond); return(NULL)})
        close(con)
        dict[[sear]] <- temp
        A2 <- temp
      }
    }
    A <- crossing(A1, A2)
    A <- t(apply(A, 1, function(x){
      tmp <- apply(data.frame(regmatches(x, gregexpr("\\d{2,}", x))), 1, as.numeric)
      x[order(order(tmp[,1], tmp[,2]))]
    }))
    A <- unique(apply(A, 1, paste, collapse=","))
    ustring <- c(ustring, paste(A, collapse="|"))
  }
  return(ustring)
}

#####################################################
#'
#' Creates a Uniformat file from a dataframe
#' 
#' 
#' Creates a .unif file containing the genetic infomation from a dataframe 
#' 
#' @aliases to.uniformat
#' @usage toUniformat(typing, fname=NULL, loci=c('A', 'B', 'C', 'DRB1', 'DQB1', 'DPB1'), dict=NULL)
#'  
#' 
#' 
#' @param typing a dataframe with first column [IND] the identifiers and the allele information in 
#' paired locus columns (Locus_1, Locus_2); unknown alleles should be [NA]
#' @param fname an optional filename to save the output; if [NULL]the output will be returned instead
#' @param loci the loci to be exported in the Uniformat file
#' @param dict the list dictionary containing the correspondence of ambiguous alleles
#' 
#' @return res Uniformat formatted data
#' 
#'  
#' 
#' @author William Lemieux \email{william.lemieux@@hema-quebec.qc.ca}
#' 
#'
#'@export
###########################################
toUniformat <- function(typing, fname=NULL, loci=c('A', 'B', 'C', 'DRB1', 'DQB1', 'DPB1'), dict=NULL){
  res <- t(apply(typing, 1, pasteUnif, loci=loci, dict=dict))
  if(is.null(fname)){
    return(res)
  } else {
    con=file(fname, "wb")
    write.table(res, con, sep="\t", col.names=F, row.names=F, quote = F)
    close(con)
  }
}