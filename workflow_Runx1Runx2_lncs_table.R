library(seqsetvis)
library(rtracklayer)
library(magrittr)
library(data.table)

###load reference
source("../jrb_R_scripts/parse_gtf.dt.R")

ref_extension = function(ext,
                         np_files,
                         bw_files,
                         goi,
                         ref_dt){
    ref_dt = copy(ref_dt)
    ref_dt[strand == "+", end := start]
    ref_dt[strand == "-", start := end]
    ref_dt[,start := start - ext]
    ref_dt[,end := end + ext - 1]
    ref_gr = GRanges(ref_dt)

    #filter to desired lists here
    ref_gr = subset(ref_gr, gene_id %in% goi)

    ### load peak calls

    #from: https://charlesjb.github.io/How_to_import_narrowPeak/
    extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                              qValue = "numeric", peak = "integer")
    np_grs <- lapply(np_files, function(f){
        import(f, format = "BED")
    })

    ### perform 3-way intersection
    # np_overlap = overlapIntervalSets(np_grs)
    np_overlap = overlapIntervalSets(append(list("lnc" = ref_gr), np_grs), use_first = T)
    np_overlap$no_hit = NULL
    names(np_overlap) = ref_gr$gene_id
    return(np_overlap)
}

np_files = c("/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_Runx1_pooled/MDA231_Runx1_pooled_peaks_passIDR.05_FEgt6.bed",
             "/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_Runx2_pooled/MDA231_Runx2_pooled_peaks_passIDR.05_FEgt6.bed")
names(np_files) = np_files %>%
    basename() %>%
    sub("GSE98551_", "", x = .) %>%
    sub("_pooled_peaks_passIDR.05_FEgt6.bed", "", x = .)

bw_files = c("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MDA231_MK_runx/MDA231_Runx1_pooled_FE.bw",
             "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MDA231_MK_runx/MDA231_Runx2_pooled_FE.bw")

# bw_files = https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE98551&format=file&file=GSE98551%5FMCF10AT1%5FCTCF%5Fpooled%5FFE%2Ebw
names(bw_files) = bw_files %>%
    basename() %>%
    sub("GSE98551_", "", x = .) %>%
    sub("_pooled_FE.bw", "", x = .)

goi = c("ENSG00000281202.2", "ENSG00000259932.1", "ENSG00000238012.1",
        "ENSG00000250626.2", "ENSG00000224137.1", "ENSG00000251359.4",
        "ENSG00000226520.1", "ENSG00000215483.9", "ENSG00000264954.1",
        "ENSG00000261114.1", "ENSG00000260737.5", "ENSG00000243485.4",
        "ENSG00000279245.1", "ENSG00000272681.2")

ref_dt = parse_gtf.dt("/slipstream/home/joeboyd/gencode.v25.annotation.gtf", additional_attribs = "gene_type")
ref_dt = ref_dt[gene_type %in% c("lincRNA", "antisense_RNA", "sense_intronic") & seqnames != "chrY"]

olap_10kb = ref_extension(10*10^3, np_files = np_files, bw_files = bw_files, goi = goi, ref_dt = ref_dt)
olap_30kb = ref_extension(30*10^3, np_files = np_files, bw_files = bw_files, goi = goi, ref_dt = ref_dt)
olap_100kb = ref_extension(100*10^3, np_files = np_files, bw_files = bw_files, goi = goi, ref_dt = ref_dt)

