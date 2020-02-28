context("EpiTxDb.Hs.hg38")
test_that("EpiTxDb.Hs.hg38:",{
    etdb <- EpiTxDb.Hs.hg38.RMBase()
    expect_s4_class(etdb,"EpiTxDb")
})
