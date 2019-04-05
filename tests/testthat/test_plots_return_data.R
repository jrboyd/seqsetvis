testthat::context("DataFromPlots")
#test return_data functionality for plotting functions.  added 10/23/18, v1.1.6

library(seqsetvis)
library(testthat)
library(GenomicRanges)

test_signals = CTCF_in_10a_profiles_dt
test_features = CTCF_in_10a_overlaps_gr


for (ssvFeatureFUN_name in c("ssvFeatureBars", "ssvFeatureEuler", "ssvFeaturePie", "ssvFeatureBinaryHeatmap", "ssvFeatureVenn")){
    test_that(paste(ssvFeatureFUN_name, "return_data"), {
        ssvFeatureFUN = get(ssvFeatureFUN_name)
        p_default = ssvFeatureFUN(test_features)
        p_dataFALSE = ssvFeatureFUN(test_features, return_data = FALSE)
        p_dataTRUE = ssvFeatureFUN(test_features, return_data = TRUE)

        expect_s3_class(p_default, "ggplot")
        expect_s3_class(p_dataFALSE, "ggplot")
        expect_s3_class(p_dataTRUE, "data.table")
    })

}

for (ssvSignalFUN_name in c("ssvSignalBandedQuantiles", "ssvSignalHeatmap", "ssvSignalLineplot", "ssvSignalLineplotAgg")){
    test_that(paste(ssvFeatureFUN_name, "return_data"), {
        ssvSignalFUN = get(ssvSignalFUN_name)
        p_default = ssvSignalFUN(test_signals)
        p_dataFALSE = ssvSignalFUN(test_signals, return_data = FALSE)
        p_dataTRUE = ssvSignalFUN(test_signals, return_data = TRUE)

        expect_s3_class(p_default, "ggplot")
        expect_s3_class(p_dataFALSE, "ggplot")
        expect_s3_class(p_dataTRUE, "data.table")
    })
}

ssvSignalFUN_name = "ssvSignalScatterplot"
test_that(paste(ssvFeatureFUN_name, "return_data"), {
    ssvSignalFUN = get(ssvSignalFUN_name)
    p_default = ssvSignalFUN(test_signals, x_name = "MCF10A_CTCF", y_name = "MCF10AT1_CTCF")
    p_dataFALSE = ssvSignalFUN(test_signals, x_name = "MCF10A_CTCF", y_name = "MCF10AT1_CTCF", return_data = FALSE)
    p_dataTRUE = ssvSignalFUN(test_signals, x_name = "MCF10A_CTCF", y_name = "MCF10AT1_CTCF", return_data = TRUE)

    expect_s3_class(p_default, "ggplot")
    expect_s3_class(p_dataFALSE, "ggplot")
    expect_s3_class(p_dataTRUE, "data.table")
})
