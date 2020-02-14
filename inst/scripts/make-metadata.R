
# Base data for all data sets --------------------------------------------------

df_Base <- DataFrame(
  BiocVersion = "3.11",
  SourceVersion = NA,
  Coordinate_1_based = TRUE,
  Maintainer = "Felix G.M. Ernst <felix.gm.ernst@outlook.com>"
)

RMBaseURL <- "http://rna.sysu.edu.cn/rmbase/"
snoRNAdbURL <- "https://www-snorna.biotoul.fr/"
tRNAdbURL <- "http://trna.bioinf.uni-leipzig.de/"

df <- rbind(
  cbind(df_Base,
        DataFrame(Title = "EpiTxDb RMBase v2.0 for Homo sapiens hg38", 
                  Description = paste0(
                    ""), 
                  SourceType = "TXT",
                  SourceUrl = RMBaseURL,
                  DataProvider = "RMBase v2.0",
                  RDataClass = "SQLiteFile", 
                  DispatchClass = "SQLiteFile",
                  RDataPath = "EpiTxDb.Hs.hg38/EpiTxDb.Hs.hg38.RMBase.sqlite")),
  cbind(df_Base,
        DataFrame(Title = "EpiTxDb snoRNAdb for Homo sapiens hg38", 
                  Description = paste0(
                    ""),
                  SourceType = "CSV",
                  SourceUrl = snoRNAdbURL,
                  DataProvider = "snoRNAdb",
                  RDataClass = "SQLiteFile", 
                  DispatchClass = "SQLiteFile",
                  RDataPath = "EpiTxDb.Hs.hg38/EpiTxDb.Hs.hg38.snoRNAdb.sqlite")),
  cbind(df_Base,
        DataFrame(Title = "EpiTxDb tRNAdb for Homo sapiens hg38", 
                  Description = paste0(
                    ""),
                  SourceType = "TXT",
                  SourceUrl = tRNAdbURL,
                  DataProvider = "tRNAdb",
                  RDataClass = "SQLiteFile", 
                  DispatchClass = "SQLiteFile",
                  RDataPath = "EpiTxDb.Hs.hg38/EpiTxDb.Hs.hg38.tRNAdb.sqlite"))
)

df$Species <- "Homo sapiens"
df$TaxonomyId <- "9606"
df$SourceVersion <- Sys.time()
df$Genome <- "hg38"
df$Tags <- "EpiTxDb:hg38:Modification:Epitranscriptomics"

write.csv(df, file = "inst/extdata/metadata.csv", row.names = FALSE)
