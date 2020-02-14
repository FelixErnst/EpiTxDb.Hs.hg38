#' @title Annotation package for \code{EpiTxDb} objects
#'
#' @author Felix G M Ernst [aut]
#'
#' @description
#' This package loads one or more \code{EpiTxDb} objects. Such \code{EpiTxDb}
#' objects are an R interface to prefabricated databases contained by this
#' package.
#' 
#' The names of any objects exposed by this package indicate the origin and
#' resources exposed. So for example \code{EpiTxDb.Hs.hg38.snoRNAdb} would be a
#' \code{EpiTxDb} object of Homo sapiens data from snoRNAdb build based on the
#' hg38 genome build.
#'
#' @return a \code{\link[EpiTxDb:EpiTxDb-class]{EpiTxDb}} object 
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
#' @name EpiTxDb.Hs.hg38
#' 
#' @examples 
#' EpiTxDb.Hs.hg38.snoRNAdb
NULL

#' @import EpiTxDb
NULL

.check_version <- function(version){
  if(!is.character(version) || length(version) != 1L){
    stop("'version' must be single character value.",call. = FALSE)
  }
  if(!(version %in% AH_DATA$version)){
    stop("'version' must be valid version. Currently valid versions are: '",
         paste(unique(AH_DATA$version), collapse = "', '"),"'",
         call. = FALSE)
  }
}

#' @rdname EpiTxDb.Hs.hg38
#' @export
EpiTxDb.Hs.hg38.RMBase <- function(version = "1"){
  .check_version(version)
  ah <- AnnotationHub()
  resource <- ah[[AH_DATA[AH_DATA$version == version,"RMBase"]]]
  return(resource)
}

#' @rdname EpiTxDb.Hs.hg38
#' @export
EpiTxDb.Hs.hg38.snoRNAdb <- function(version = "1"){
  .check_version(version)
  ah <- AnnotationHub()
  resource <- ah[[AH_DATA[AH_DATA$version == version,"snoRNAdb"]]]
  return(resource)
}

#' @rdname EpiTxDb.Hs.hg38
#' @export
EpiTxDb.Hs.hg38.tRNAdb <- function(version = "1"){
  .check_version(version)
  ah <- AnnotationHub()
  resource <- ah[[AH_DATA[AH_DATA$version == version,"tRNAdb"]]]
  return(resource)
}

# version information ----------------------------------------------------------

AH_DATA <- data.frame(version = "1",
                      RMBase = "AH00000",
                      snoRNAdb = "AH00000",
                      tRNAdb = "AH00000")

# AH_DATA <- rbind(AH_DATA,
#                  data.frame(version = "1.0",
#                             RMBase = "AH00000",
#                             snoRNAdb = "AH00000",
#                             tRNAdb = "AH00000"))
