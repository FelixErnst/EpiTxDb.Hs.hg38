context("EpiTxDb.Hs.hg38")
test_that("EpiTxDb.Hs.hg38:",{
    path <- system.file("extdata", package = "EpiTxDb.Hs.hg38")
    files <- dir(path)
    files <- files[grepl("*\\.sqlite", files)]
    expect_equal(length(files), 3L)
})
