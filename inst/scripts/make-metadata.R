
# Base data for all data sets --------------------------------------------------
library(S4Vectors)

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
                    "Information from the RMBase v2.0 was downloaded and ",
                    "imported as EpiTxDb database. All valid modification ",
                    "types for Homo sapiens/hg19 were used, lift over to hg38,",
                    " checked the nucleotide sequence."), 
                  SourceType = "BED",
                  SourceUrl = RMBaseURL,
                  DataProvider = "RMBase v2.0",
                  RDataClass = "SQLiteFile", 
                  DispatchClass = "SQLiteFile",
                  RDataPath = "EpiTxDb.Hs.hg38/EpiTxDb.Hs.hg38.RMBase.sqlite")),
  cbind(df_Base,
        DataFrame(Title = "EpiTxDb snoRNAdb for Homo sapiens hg38", 
                  Description = paste0(
                    "Information from the snoRNAdb was downloaded, manually ",
                    "adjusted for changes in recent rRNA sequences from ",
                    "H. sapiens and stored as EpiTxDb database. The ",
                    "information provided match hg38 release sequences."),
                  SourceType = "CSV",
                  SourceUrl = snoRNAdbURL,
                  DataProvider = "snoRNAdb",
                  RDataClass = "SQLiteFile", 
                  DispatchClass = "SQLiteFile",
                  RDataPath = "EpiTxDb.Hs.hg38/EpiTxDb.Hs.hg38.snoRNAdb.sqlite")),
  cbind(df_Base,
        DataFrame(Title = "EpiTxDb tRNAdb for Homo sapiens hg38", 
                  Description = paste0(
                    "tRNAdb entries for RNA database were imported using ",
                    "tRNAdbImport and the modification information extracted. ",
                    "tRNA sequences were matched against transcripts, ",
                    "results combined and stored as EpiTxDb database."),
                  SourceType = "XML",
                  SourceUrl = tRNAdbURL,
                  DataProvider = "tRNAdb",
                  RDataClass = "SQLiteFile", 
                  DispatchClass = "SQLiteFile",
                  RDataPath = "EpiTxDb.Hs.hg38/EpiTxDb.Hs.hg38.tRNAdb.sqlite")),
  cbind(df_Base,
        DataFrame(Title = "Chain file for Homo sapiens rRNA hg19 to hg38", 
                  Description = paste0(
                    ""),
                  SourceType = "Chain",
                  SourceUrl = "https://www.ncbi.nlm.nih.gov/nuccore/U13369;https://www.ncbi.nlm.nih.gov/nuccore/NR_003287.4;https://www.ncbi.nlm.nih.gov/nuccore/NR_003286.4;https://www.ncbi.nlm.nih.gov/nuccore/NR_003285.3",
                  DataProvider = "NCBI",
                  RDataClass = "ChainFile", 
                  DispatchClass = "ChainFile",
                  RDataPath = "EpiTxDb.Hs.hg38/rRNA.hg19Tohg38.liftOver")),
  cbind(df_Base,
        DataFrame(Title = "Chain file for Homo sapiens rRNA hg38 to hg19", 
                  Description = paste0(
                    ""),
                  SourceType = "Chain",
                  SourceUrl = "https://www.ncbi.nlm.nih.gov/nuccore/U13369;https://www.ncbi.nlm.nih.gov/nuccore/NR_003287.4;https://www.ncbi.nlm.nih.gov/nuccore/NR_003286.4;https://www.ncbi.nlm.nih.gov/nuccore/NR_003285.3",
                  DataProvider = "NCBI",
                  RDataClass = "ChainFile", 
                  DispatchClass = "ChainFile",
                  RDataPath = "EpiTxDb.Hs.hg38/rRNA.hg38Tohg19.liftOver"))
)

df$Species <- "Homo sapiens"
df$TaxonomyId <- "9606"
df$SourceVersion <- Sys.time()
df$Genome <- "hg38"
df$Tags <- "EpiTxDb:hg38:Modification:Epitranscriptomics"

write.csv(df, file = "inst/extdata/metadata.csv", row.names = FALSE)
