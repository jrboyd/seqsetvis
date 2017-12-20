library(peakvisr)
library(testthat)
set_a = 1:8
set_b = 5:9
set_c = 8:10

sets_list = list(set_a, set_b, set_c)
sets_list_named = list("a" = set_a, "b" = set_b, "c" = set_c)

memb_table = set_list2memb(sets_list)

gr_a = GRanges("chr1", IRanges(set_a*10, set_a*10+5))
gr_b = GRanges("chr1", IRanges(set_b*10-2, set_b*10+3))
gr_c = GRanges("chr1", IRanges(set_c*10+-4, set_c*10+1))

sets_gr_list = list(gr_a, gr_b, gr_c)
sets_gr_list_named = list("a" = gr_a, "b" = gr_b, "c" = gr_c)

gr_olap = overlapIntervalSets(sets_gr_list)
gr_olap_named = overlapIntervalSets(sets_gr_list_named)

setPlotMakeMT(sets_list)
setPlotMakeMT(sets_list_named)

setPlotMakeMT(sets_gr_list)
setPlotMakeMT(sets_gr_list_named)

setPlotMakeMT(as.matrix(memb_table))
setPlotMakeMT(memb_table)

setPlotMakeMT(gr_olap)
setPlotMakeMT(mcols(gr_olap))
setPlotMakeMT(as.matrix(mcols(gr_olap)))
setPlotMakeMT(as.data.frame(mcols(gr_olap)))
setPlotMakeMT(as.data.table(mcols(gr_olap)))

setPlotMakeMT(gr_olap_named)
setPlotMakeMT(mcols(gr_olap_named))
setPlotMakeMT(as.matrix(mcols(gr_olap_named)))
setPlotMakeMT(as.data.frame(mcols(gr_olap_named)))
setPlotMakeMT(as.data.table(mcols(gr_olap_named)))

setPlotMakeMT(GRangesList(sets_gr_list))
