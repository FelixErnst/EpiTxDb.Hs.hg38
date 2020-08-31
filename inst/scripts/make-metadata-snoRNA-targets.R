
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
        DataFrame(Title = "Sequences of snoRNA targets of Homo sapiens hg38", 
                  Description = paste0(
                    "Fasta file for snoRNA targets based on genomic sequences ",
                    "for Homo sapiens hg38."),
                  SourceType = "FASTA",
                  SourceUrl = "https://www.ncbi.nlm.nih.gov/gene",
                  DataProvider = "NCBI",
                  RDataClass = "FaFile", 
                  DispatchClass = "FaFile",
                  RDataPath = "EpiTxDb.Hs.hg38/snoRNA.targets.hg38.fa"))
)

df$Species <- "Homo sapiens"
df$TaxonomyId <- "9606"
df$SourceVersion <- Sys.time()
df$Genome <- "hg38"
df$Tags <- "EpiTxDb:hg38:Modification:Epitranscriptomics"

write.csv(df, file = "inst/extdata/metadata-snoRNA-targets.csv",
          row.names = FALSE)
