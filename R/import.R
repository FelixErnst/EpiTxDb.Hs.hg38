
import_from_snoRNAdb <- function(snoRNAdb, orgdb){
  mod_id <- seq_len(nrow(snoRNAdb))
  mod_name <- paste0(snoRNAdb$modification,"_",snoRNAdb$position)
  mod_type <- snoRNAdb$modification
  mod_start <- snoRNAdb$position
  mod_end <- snoRNAdb$position
  modifications <- data.frame(mod_id = mod_id,
                              mod_name = mod_name,
                              mod_type = mod_type,
                              mod_start = mod_start,
                              mod_end = mod_end,
                              stringsAsFactors = FALSE)
  
  gene <- select(orgdb,keys = "fibrillarin",
                 columns = c("GENENAME","ENSEMBL","ENTREZID","ENZYME"),
                 keytype = "GENENAME")
  gene <- gene[gene$ENSEMBL == "ENSG00000105202",]
  mod_rank <- 1L
  mod_type <- snoRNAdb$modification
  genename <- gene$GENENAME
  ensembl <- gene$ENSEMBL
  ensembltrans <- ""
  entrezid <- gene$ENTREZID
  enzyme <- gene$ENZYME
  reactions <- data.frame(mod_id = mod_id,
                          mod_rank = mod_rank,
                          mod_type = mod_type,
                          reaction_genename = genename,
                          reaction_ensembl = ensembl,
                          reaction_ensembltrans = ensembltrans,
                          reaction_entrezid = entrezid,
                          reaction_enzyme = enzyme,
                          stringsAsFactors = FALSE)
  
  specifier_genename <- snoRNAdb$guide
  specifier_f <- specifier_genename != "unknown"
  specifier_genename <- strsplit(as.character(specifier_genename),",")[specifier_f]
  specifier_lengths <- lengths(specifier_genename)
  specifier_type <- "snoRNA"
  specifier_mod_id <- unlist(Map(rep,mod_id[specifier_f],specifier_lengths))
  specifier_entrezid <- mapIds(orgdb,unlist(specifier_genename),"ENTREZID",
                               "SYMBOL", multiVals = "list")
  specifier_ensembl <- mapIds(orgdb,unlist(specifier_genename),"ENSEMBL",
                               "SYMBOL", multiVals = "list")
  specifier_ensembl_lengths <- lengths(specifier_ensembl)
  specifier_mod_id <- unlist(Map(rep,specifier_mod_id,specifier_ensembl_lengths))
  specifier_genename <- unlist(Map(rep,unlist(specifier_genename),specifier_ensembl_lengths))
  specifier_entrezid <- unlist(Map(rep,specifier_entrezid,specifier_ensembl_lengths))
  
  specifier <- data.frame(mod_id = specifier_mod_id,
                          specifier_type = specifier_type,
                          specifier_genename = unlist(specifier_genename),
                          specifier_entrezid = specifier_entrezid,
                          specifier_ensembl = specifier_orgdb$ENSEMBL,
                          stringsAsFactors = FALSE)
  
  transcripts_entrez_id <- select(orgdb,as.character(snoRNAdb$hgnc_symbol),
                                  c("ENTREZID","SYMBOL"),"SYMBOL")[,"ENTREZID"]
  transcripts <- data.frame(mod_id = mod_id,
                            entrezid = transcripts_entrez_id,
                            stringsAsFactors = FALSE)
  makeTxModDb(modifications, reactions, specifier, transcripts)
}
