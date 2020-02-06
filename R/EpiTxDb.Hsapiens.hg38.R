#' @title Annotation package for EpiTxDb object(s)
#'
#' @author Felix G M Ernst [aut]
#'
#' @description
#' This package loads one or more EpiTxDb objects. Such EpiTxDb objects are an
#' R interface to prefabricated databases contained by this package.
#' 
#' The names of any objects exposed by this package indicate the origin and
#' resources exposed.  So for example EpiTxDb.Hsapiens.hg38.snoRNAdb would be a
#' EpiTxDb object of Homo sapiens data from snoRNAdb build based on the hg38
#' build.
#'
#' @note 
#' 
#' @seealso
#' \itemize{
#' \item{\code{\link[EpiTxDb:modifications]{modifications}}}
#' \item{\code{\link[EpiTxDb:modifications]{reactions}}}
#' \item{\code{\link[EpiTxDb:modifications]{specifies}}}
#' \item{\code{\link[EpiTxDb:modifications]{modificationsByTranscript}}}
#' \item{\code{\link[EpiTxDb:modifications]{modifiedSeqsByTranscript}}}
#' }
#' 
#' @docType package
#' @name EpiTxDb.Hsapiens.hg38
#' 
#' @examples 
#' library(EpiTxDb.Hsapiens.hg38)
#' EpiTxDb.Hsapiens.hg38.snoRNAdb
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
