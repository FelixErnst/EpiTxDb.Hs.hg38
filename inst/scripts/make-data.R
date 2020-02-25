# get annotation data
library(AnnotationHub)
library(BSgenome.Hsapiens.UCSC.hg38)
library(org.Hs.eg.db)
library(EpiTxDb)
library(RSQLite)
library(GenomicRanges)

# get annotation hub
ah <- AnnotationHub()

# get ensembl annotation
edb <- query(ah, c("EnsDb","Homo sapiens","99"))[[1]]
seqlevelsStyle(edb) <- "UCSC"
seqlevels <- paste0("chr",c(seq_len(22),"Y","X","M"))

# get hg19 to hg38 chain file
chainfiles <- query(ah , c("hg38", "hg19", "chainfile"))
chain <- chainfiles[[2]]

# get transcript annotations from ensemble db
assemble_tx <- function(edb, genome){
  .split_tRNA_by_intron <- function(gr){
    if(length(gr) > 1L){
      stop(".")
    }
    regres <- regmatches(gr$intron,regexec(".*pos ([0-9]+)-([0-9]+)<BR>",gr$intron))
    start <- as.integer(regres[[1L]][2L])
    end <- as.integer(regres[[1L]][3L])
    if(as.logical(strand(gr) == "-")){
      ranges <- IRanges::IRanges(start = c(end(gr)-start+2L,start(gr)),
                                 end = c(end(gr),end(gr)-end))
    } else {
      ranges <- IRanges::IRanges(start = c(start(gr),start(gr)+end),
                                 end = c(start(gr)+start-2L,end(gr)))
    }
    ans <- GenomicRanges::GRanges(seqnames = seqnames(gr),
                                  ranges = ranges,
                                  strand = strand(gr),
                                  mcols(gr))
    mcols(ans) <- mcols(ans)[,colnames(mcols(ans)) != "intron",drop=FALSE]
    ans
  }
  # get exons
  tx <- exonsBy(edb,"tx")
  tx_id <- IRanges::CharacterList(Map(rep,names(tx),lengths(tx)))
  mcols(tx, level="within")[,"tx_id"] <- tx_id
  genome(tx) <- genome
  
  # get tRNA annotation
  FDb.tRNAs <- makeFeatureDbFromUCSC(genome,"tRNAs","tRNAs")
  tRNAs <- features(FDb.tRNAs)
  mcols(tRNAs) <- mcols(tRNAs)[,c("intron"),drop=FALSE]
  mcols(tRNAs)$tx_id <- names(tRNAs)
  names(tRNAs) <- NULL
  tRNAs <- split(tRNAs,mcols(tRNAs)$tx_id)
  # incorporate intron annotations
  has_intron <- mcols(tRNAs,level="within")[,"intron"] != "No tRNA introns"
  has_intron <- unlist(unique(has_intron))
  tRNAs[has_intron] <- GenomicRanges::GRangesList(lapply(tRNAs[has_intron],.split_tRNA_by_intron))
  exon_id <- relist(paste(unlist(Map(rep,names(tRNAs),lengths(tRNAs))),
                          unlist(lapply(lengths(tRNAs),seq_len)),
                          sep="_"),
                    PartitioningByEnd(tRNAs))
  mcols(tRNAs, level="within")[,"exon_id"] <- exon_id
  #
  mcols(tx, level="within") <- mcols(tx, level="within")[,c("tx_id","exon_id")]
  mcols(tRNAs, level="within") <- mcols(tRNAs, level="within")[,c("tx_id","exon_id")]
  ans <- c(tx,tRNAs)
  # fix the seqinfo
  seqlevels <- paste0("chr",c(seq_len(22),"Y","X","M"))
  ans <- ans[seqnames(ans) %in% seqlevels]
  seqlevels(ans) <- seqlevels
  circ <- isCircular(seqinfo(ans))
  circ[names(circ) == "chrM"] <- TRUE
  isCircular(seqinfo(ans)) <- circ
  #
  ans[lengths(ans) != 0L]
}

tx <- assemble_tx(edb, "hg38")

################################################################################
# functions for import
################################################################################

import.RMBase <- function(bs, organism, genome, type, chain){
  metadata <- data.frame(name = c("Data source","Organism","Genome",
                                  "Coordinates"),
                         value = c("RMBase v2.0","Homo sapiens","hg38",
                                   "per Genome"))
  #
  files <- downloadRMBaseFiles(organism, genome, type)
  gr <- getRMBaseDataAsGRanges(files)
  # use liftOver to get the hg38 coordinates
  gr <- unlist(liftOver(gr,chain))
  f <- !duplicated(paste0(as.character(gr),"-",
                          mcols(gr)$mod_type))
  gr <- gr[f]
  seq <- getSeq(bs, seqlevels(gr))
  seq_rna <- as(seq, "RNAStringSet")
  colnames(mcols(gr)) <- gsub("mod_type","mod",colnames(mcols(gr)))
  # do plus and minus strand separatly, since removeIncompatibleModifications
  # only accepts plus strand
  gr_plus <- 
    Modstrings::removeIncompatibleModifications(gr[strand(gr) == "+"], seq_rna)
  gr_minus <- gr[strand(gr) == "-"]
  strand(gr_minus) <- "+"
  gr_minus <-
    Modstrings::removeIncompatibleModifications(gr_minus, complement(seq_rna))
  strand(gr_minus) <- "-"
  gr <- c(gr_plus,gr_minus)
  gr <- gr[order(seqnames(gr),start(gr),strand(gr))]
  colnames(mcols(gr)) <- gsub("^mod$","mod_type",colnames(mcols(gr)))
  metadata <- EpiTxDb:::.add_sequence_check_to_metadata(metadata)
  #
  mcols(gr)$mod_id <- seq_along(gr)
  makeEpiTxDbFromGRanges(gr, metadata = metadata)
}

import_from_tRNAdb <- function(organism, bs, tx){
  metadata <- data.frame(name = c("Data source","Organism","Genome",
                                  "Coordinates"),
                         value = c("tRNAdb","Homo sapiens","hg38",
                                   "per Transcript"))
  #
  seq <- getSeq(bs,tx)
  seq <- relist(unlist(unlist(seq)),
                IRanges::PartitioningByWidth(sum(nchar(seq))))
  seq_rna <- as(seq,"RNAStringSet")
  gr <- gettRNAdbDataAsGRanges(organism, sequences = seq_rna)
  gr <- gr[!duplicated(paste0(as.character(gr),"-",gr$mod_type))]
  # fix an error in the tRNAdb for tRNA Tyr at position 20
  fix <- split(gr[start(gr) == 20],seqnames(gr[start(gr) == 20]))
  fix_f <- lengths(fix) == 2L
  fix_replace <- fix[fix_f][mcols(fix[fix_f], level="within")[,"mod_type"] == "acp3U"]
  fix[fix_f] <- fix_replace
  gr <- c(gr[start(gr) != 20],unlist(fix))
  gr <- gr[order(gr)]
  names(gr) <- NULL
  #
  colnames(mcols(gr)) <- gsub("mod_type","mod",colnames(mcols(gr)))
  gr <- Modstrings::removeIncompatibleModifications(gr, seq_rna)
  colnames(mcols(gr)) <- gsub("^mod$","mod_type",colnames(mcols(gr)))
  metadata <- EpiTxDb:::.add_sequence_check_to_metadata(metadata)
  #
  makeEpiTxDbFromGRanges(gr, metadata = metadata)
}

import_from_snoRNAdb <- function(snoRNAdb, orgdb){
  metadata <- data.frame(name = c("Data source","Organism","Genome",
                                  "Coordinates"),
                         value = c("snoRNAdb","Homo sapiens","hg38",
                                   "per Transcript"))
  # Modifications
  mod_id <- seq_len(nrow(snoRNAdb))
  mod_name <- paste0(snoRNAdb$modification,"_",snoRNAdb$position)
  mod_type <- snoRNAdb$modification
  mod_start <- snoRNAdb$position
  mod_end <- snoRNAdb$position
  
  transcripts <- select(orgdb,as.character(snoRNAdb$hgnc_symbol),
                        c("REFSEQ","SYMBOL","ENTREZID"),"SYMBOL")
  
  modifications <- data.frame(mod_id = mod_id,
                              mod_name = mod_name,
                              mod_type = mod_type,
                              mod_start = mod_start,
                              mod_end = mod_end,
                              sn_id = as.integer(transcripts$ENTREZID),
                              sn_name = transcripts$REFSEQ,
                              stringsAsFactors = FALSE)
  
  # Reactions
  gene_fbl <- select(orgdb,keys = "FBL",
                     columns = c("GENENAME","ENSEMBL","ENTREZID"),
                     keytype = "SYMBOL")
  gene_fbl <- gene_fbl[gene_fbl$ENSEMBL == "ENSG00000105202",]
  
  gene_dkc <- select(orgdb,keys = "DKC1",
                     columns = c("GENENAME","ENSEMBL","ENTREZID"),
                     keytype = "SYMBOL")
  
  rx_rank <- 1L
  mod_type <- snoRNAdb$modification
  genename <- character(length(mod_type))
  ensembl <- character(length(mod_type))
  ensembltrans <- character(length(mod_type))
  entrezid <- character(length(mod_type))
  
  genename[mod_type == "Y"] <- gene_dkc$GENENAME
  genename[mod_type != "Y"] <- gene_fbl$GENENAME
  
  ensembl[mod_type == "Y"] <- gene_dkc$ENSEMBL
  ensembl[mod_type != "Y"] <- gene_fbl$ENSEMBL
  
  entrezid[mod_type == "Y"] <- gene_dkc$ENTREZID
  entrezid[mod_type != "Y"] <- gene_fbl$ENTREZID
  
  reactions <- data.frame(mod_id = mod_id,
                          rx_genename = genename,
                          rx_rank = rx_rank,
                          rx_ensembl = ensembl,
                          rx_ensembltrans = ensembltrans,
                          rx_entrezid = entrezid,
                          stringsAsFactors = FALSE)
  
  # Specifiers
  specifier_genename <- snoRNAdb$guide
  specifier_f <- specifier_genename != "unknown"
  specifier_genename <- strsplit(as.character(specifier_genename),",")[specifier_f]
  specifier_lengths <- lengths(specifier_genename)
  specifier_type <- "snoRNA"
  specifier_mod_id <- unlist(Map(rep,mod_id[specifier_f],specifier_lengths))
  specifier_entrezid <- mapIds(orgdb,unlist(specifier_genename),"ENTREZID",
                               "SYMBOL")
  specifier_ensembl <- mapIds(orgdb,unlist(specifier_genename),"ENSEMBL",
                              "SYMBOL")
  
  specifiers <- data.frame(mod_id = specifier_mod_id,
                           spec_type = specifier_type,
                           spec_genename = unlist(specifier_genename),
                           spec_ensembl = specifier_ensembl,
                           spec_entrezid = specifier_entrezid,
                           stringsAsFactors = FALSE)
  # References
  references <- data.frame(mod_id = mod_id,
                           ref_type = "PMID",
                           ref = "16381836")
  makeEpiTxDb(modifications, reactions, specifiers, references,
              metadata = metadata)
}

# start the import RMBase, snoRNAdb and tRNAdb data
start.import <- function(bs, orgdb, tx, chain){
  etdb <- import_from_snoRNAdb(read.csv2(RNAmodR.Data::RNAmodR.Data.snoRNAdb()),
                               orgdb)
  db <- dbConnect(SQLite(), "hub/EpiTxDb.Hs.hg38.snoRNAdb.sqlite")
  sqliteCopyDatabase(etdb$conn, db)
  dbDisconnect(etdb$conn)
  dbDisconnect(db)
  
  etdb <- import_from_tRNAdb("Homo sapiens", bs, tx)
  db <- dbConnect(SQLite(), "hub/EpiTxDb.Hs.hg38.tRNAdb.sqlite")
  sqliteCopyDatabase(etdb$conn, db)
  dbDisconnect(etdb$conn)
  dbDisconnect(db)
  
  etdb <- import.RMBase(bs, "human", "hg19",
                        listAvailableModFromRMBase("human", "hg19"),
                        chain)
  db <- dbConnect(SQLite(), "hub/EpiTxDb.Hs.hg38.RMBase.sqlite")
  sqliteCopyDatabase(etdb$conn, db)
  dbDisconnect(etdb$conn)
  dbDisconnect(db)
  return(TRUE)
}

start.import(BSgenome.Hsapiens.UCSC.hg38, org.Hs.eg.db, tx, chain)
