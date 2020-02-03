#' @title Annotation package for TxModDb object(s)
#'
#' @author Felix G M Ernst [aut]
#'
#' @description
#' This package loads one or more TxModDb objects. Such TxModDb objects are an
#' R interface to prefabricated databases contained by this package.
#' 
#' The names of any objects exposed by this package indicate the origin and
#' resources exposed.  So for example TxModDb.Hsapiens.hg38.snoRNAdb would be a
#' TxModDb object of Homo sapiens data from snoRNAdb build based on the hg38
#' build.
#'
#' @note 
#' 
#' @seealso
#' \itemize{
#' \item{\code{\link[TxModDb:modifications]{modifications}}}
#' \item{\code{\link[TxModDb:modifications]{reactions}}}
#' \item{\code{\link[TxModDb:modifications]{specifies}}}
#' \item{\code{\link[TxModDb:modifications]{modificationsByTranscript}}}
#' \item{\code{\link[TxModDb:modifications]{modifiedSeqsByTranscript}}}
#' }
#' 
#' @docType package
#' @name TxModDb.Hsapiens.hg38
#' 
#' @examples 
#' library(TxModDb.Hsapiens.hg38)
#' TxModDb.Hsapiens.hg38.snoRNAdb
NULL

.onLoad <- function(libname, pkgname)
{
  ns <- asNamespace(pkgname)
  path <- system.file("extdata", package=pkgname, lib.loc=libname)
  files <- dir(path)
  files <- files[grepl("*\\.sqlite",files)]
  for(i in seq_len(length(files))){
    db <- loadDb(system.file("extdata", files[[i]], package=pkgname, 
                             lib.loc=libname),packageName=pkgname)
    objname <- sub(".sqlite$","",files[[i]])
    assign(objname, db, envir=ns)
    namespaceExport(ns, objname)
  }
}
