context("EpiTxDb.Sc.sacCer3")
test_that("EpiTxDb.Sc.sacCer3:",{
    path <- system.file("extdata", package = "EpiTxDb.Sc.sacCer3")
    files <- dir(path)
    files <- files[grepl("*\\.sqlite", files)]
    expect_equal(length(files), 2L)
})
