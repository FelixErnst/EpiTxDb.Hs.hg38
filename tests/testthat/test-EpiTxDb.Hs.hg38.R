context("EpiTxDb.Hs.hg38")
test_that("EpiTxDb.Hs.hg38:",{
    actual <- AnnotationHubData::makeAnnotationHubMetadata(system.file(package = "EpiTxDb.Hs.hg38"),
                                                           fileName = "../../extdata/metadata.csv")
    expect_equal(length(actual), 5L)
    
})
