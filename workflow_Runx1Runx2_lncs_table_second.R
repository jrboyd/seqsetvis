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
        if(grepl(pattern = "Peak$", f)){
            import(f, format = "BED", extraCols = extraCols_narrowPeak)
        }else{
            import(f, format = "BED")
        }

    })

    ### perform 3-way intersection
    # np_overlap = overlapIntervalSets(np_grs)
    np_overlap = overlapIntervalSets(append(list("lnc" = ref_gr), np_grs), use_first = T)
    np_overlap$no_hit = NULL
    names(np_overlap) = ref_gr$gene_id
    return(np_overlap)
}

odir = getwd()
setwd("/slipstream/galaxy/uploads/working/qc_framework")
np_files = c("output_MK_MDA231_Runx/MDA231_Runx1_pooled/MDA231_Runx1_pooled_peaks_passIDR.05_FEgt6.bed",
             "output_MK_MDA231_Runx/MDA231_Runx2_pooled/MDA231_Runx2_pooled_peaks_passIDR.05_FEgt6.bed",
             "output_AF_RUNX1_ChIP/AF-MCF10A_RUNX1_pooled/AF-MCF10A_RUNX1_pooled_peaks_passIDR.05.narrowPeak",
             "output_AF_RUNX1_ChIP/AF-MCF10AT1_RUNX1_pooled/AF-MCF10AT1_RUNX1_pooled_peaks_passIDR.05.narrowPeak",
             "output_AF_RUNX1_ChIP/AF-MCF10CA1_RUNX1_pooled/AF-MCF10CA1_RUNX1_pooled_peaks_passIDR.05.narrowPeak",
             "output_Rasim_RUNX1/MCF7_RUNX1_pooled/MCF7_RUNX1_pooled_peaks_passIDR.05.narrowPeak")
names(np_files) = np_files %>%
    basename() %>%
    sub("GSE98551_", "", x = .) %>%
    sub("_pooled_peaks_passIDR.05_FEgt6.bed", "", x = .) %>%
    sub("_pooled_peaks_passIDR.05.narrowPeak", "", x = .)

names(np_files)[1:2] = paste0("MK-", names(np_files)[1:2])
names(np_files)[6] = paste0("RB-", names(np_files)[6])


soi = c("AC013451.2", #"CTD-2207P18.2",
        "MIR503HG",
        "LINC00704",
        "TMPO-AS1",
        "LINC00963",
        "SFTA1P",
        "AC091182.1", #"RP11-150O12.1",
        "AL512274.1", #"RP11-7k24.3",
        "ZEB1-AS1")

ref_dt = parse_gtf.dt("/slipstream/home/joeboyd/gencode.v27.annotation.gtf", additional_attribs = "gene_type")
# ref_dt = ref_dt[gene_type %in% c("lincRNA", "antisense_RNA", "sense_intronic", "processed_transcript") & seqnames != "chrY"]
ref_dt[,c("stable_gene_id") := tstrsplit(gene_id, "\\.", keep = 1)]

goi = ref_dt[gene_name %in% soi]$gene_id

goi2 = c("ENSG00000170846.16
ENSG00000236830.6
         ENSG00000235123.5
         ENSG00000197308.9
         ENSG00000259527.2
         ENSG00000248663.6
         ENSG00000230487.7
         ENSG00000265702.1
         ENSG00000273270.1
         ENSG00000223808.1
         ENSG00000257337.6
         ENSG00000255310.2
         ENSG00000261662.1
         ENSG00000270720.1
         ENSG00000231028.8
         ENSG00000244055.1
         ENSG00000253106.1
         ENSG00000259673.5
         ENSG00000271971.1
         ENSG00000214770.2
         ENSG00000257167.2
         ENSG00000267751.5
         ENSG00000272524.1
         ENSG00000236256.9
         ENSG00000237523.1
         ENSG00000273132.1
         ENSG00000226808.1
         ENSG00000234477.1
         ENSG00000250159.6
         ENSG00000259070.5
         ENSG00000268621.5
         ENSG00000242147.1")
goi2 = strsplit(goi2, "\n ?")[[1]] %>% gsub(" ", "", .)
goi = c(goi, goi2)
### overlap with different extension sizes
todo = c(2, 10, 30, 100)*10^3
names(todo) = paste0(c(2, 10, 30, 100), "kb")
olaps = lapply(todo, function(ext){
    ref_extension(ext = ext, np_files = np_files, bw_files = bw_files, goi = goi, ref_dt = ref_dt)
})

setwd(odir)

### prepare results for merging
olaps_df = lapply(names(olaps), function(nam){
    olap_this = olaps[[nam]]
    df_this = mcols(olap_this)[,-1]
    colnames(df_this) = paste(colnames(df_this), nam)
    df_this$id = names(olap_this)
    df_this
})

###merge
olaps_master = olaps_df[[1]]
hidden = lapply(olaps_df[-1], function(x){
    olaps_master <<- merge(olaps_master, x, by = "id")
    NULL
})
olaps_master = as.data.table(olaps_master)
olaps_master[, c("stable_gene_id", "version_gene_id") := tstrsplit(id, "\\.")]
olaps_master$id = NULL
olaps_master$version_gene_id = NULL
o = order(rep(seq_along(np_files)-1, length(todo)))
olaps_master = olaps_master[, c(o, length(o)+1), with = F]


setwd("~/zdrive/RNA_seq_data/Breast/")
DE_files = c("MCF10A_shRunx1/DESeq2/DESeq2_EV_vs_C1shRUNX1_v25.txt",
             "MCF10A_shRunx1/DESeq2/DESeq2_EV_vs_C4shRUNX1_v25.txt",
             "MCF10A-AT1_CA1_shRUNX2/DESeq2/MCF10A-AT1_shRUNX2/DESeq2_MCF10A-AT1_EV_vs_shRUNX2_v25.txt",
             "MCF10A-AT1_CA1_shRUNX2/DESeq2/MCF10A-CA1_shRUNX2/2-way_comparison_EV_and_shRUNX2_only/DESeq2_MCF10A-CA1_EV_vs_shRUNX2_v25.txt",
             "MCF7_shRunx1/DESeq2/DESeq2_EV_vs_C1shRunx1_v25.txt",
             "MCF7_shRunx1/DESeq2/DESeq2_EV_vs_C4shRunx1_v25.txt",
             "MDA-MB-231_Runx1/DESeq2/MDA_Runx1_3replicates_siRunx1_v21/DESeq2_MDA_siNS_vs_MDA_siRunx1.txt",
             "All_5_BC_cell_Lines/v25/DESeq2_AT1_vs_CA1_v25.txt",
             "All_5_BC_cell_Lines/v25/DESeq2_AT1_vs_MCF7_v25.txt",
             "All_5_BC_cell_Lines/v25/DESeq2_AT1_vs_MDA_v25.txt",
             "All_5_BC_cell_Lines/v25/DESeq2_CA1_vs_MCF7_v25.txt",
             "All_5_BC_cell_Lines/v25/DESeq2_CA1_vs_MDA_v25.txt",
             "All_5_BC_cell_Lines/v25/DESeq2_MCF10A_vs_AT1_v25.txt",
             "All_5_BC_cell_Lines/v25/DESeq2_MCF10A_vs_CA1_v25.txt",
             "All_5_BC_cell_Lines/v25/DESeq2_MCF10A_vs_MCF7_v25.txt",
             "All_5_BC_cell_Lines/v25/DESeq2_MCF10A_vs_MDA_v25.txt",
             "All_5_BC_cell_Lines/v25/DESeq2_MCF7_vs_MDA_v25.txt")
             # "MCF10A_vs_MCF10A-AT1/DESeq2_MCF10A_vs_MCF10A-AT1_Significant.txt",
             # "MCF10A_vs_MCF10A-CA1/DESeq2_MCF10A_vs_MCF10A-CA1_Significant.txt",
             # "MCF10A-AT1_vs_MCF10A-CA1/DESeq2_MCF10A-AT1_vs_MCF10A-CA1_Significant.txt",
             # "MCF10A_MCF7_MDA231/DESeq2/DESeq2_MCF10A_vs_MCF7_Significant_DOWN_in_MCF10A.txt",
             # "MCF10A_MCF7_MDA231/DESeq2/DESeq2_MCF10A_vs_MCF7_Significant_UP_in_MCF10A.txt",
             # "MCF10A_MCF7_MDA231/DESeq2/DESeq2_MCF10A_vs_MDA-MB-231_Significant_DOWN_in_MCF10A.txt",
             # "MCF10A_MCF7_MDA231/DESeq2/DESeq2_MCF10A_vs_MDA-MB-231_Significant_UP_in_MCF10A.txt",
             # "MCF10A_MCF7_MDA231/DESeq2/DESeq2_MCF7_vs_MDA-MB-231_Significant_DOWN_in_MCF7.txt",
             # "MCF10A_MCF7_MDA231/DESeq2/DESeq2_MCF7_vs_MDA-MB-231_Significant_UP_in_MCF7.txt")
names(DE_files) = DE_files %>%
    basename() %>%
    sub("DESeq2_", "", .) %>%
    sub("_v25.txt", "", .) %>%
    sub(".txt", "", .) %>%
    sub("_Significant", "", .)
names(DE_files)[c(1:2, 5:6)] = paste0(c("MCF10A_", "MCF10A_", "MCF7_", "MCF7_"), names(DE_files)[c(1:2, 5:6)], sep = "_")
names(DE_files) = sub("_$", "", names(DE_files))
names(DE_files)[3] = paste0(names(DE_files)[3], "shRUNX2")

DE_df = lapply(DE_files, function(x)read.table(x, row.names = 1, header = T))

setwd(odir)

DE_dt = lapply(DE_df, function(x)as.data.table(x, keep.rownames= T))
DE_dt = rbindlist(DE_dt, idcol = "comparison")
DE_dt[, c("stable_gene_id", "version_gene_id") := tstrsplit(rn, "\\.")]
DE_dt[, .N, by = comparison]
cnts = DE_dt[, .N, by = stable_gene_id]
hist(cnts$N)

DE_dt$comparison = factor(DE_dt$comparison, levels = names(DE_files))

fc_df = as.data.frame(dcast(DE_dt, stable_gene_id ~ comparison, value.var = "log2FoldChange"))
padj_dt = dcast(DE_dt, stable_gene_id ~ comparison, value.var = "padj")

is_sig = padj_dt[,-1] < .05
fc_df[,-1][!as.data.frame(is_sig)] = NA

cn = colnames(fc_df)[-1]
tos = sapply(cn, function(x)strsplit(x, "_vs_")[[1]][1])
froms = sapply(cn, function(x)strsplit(x, "_vs_")[[1]][2])
cn = paste("from",froms, "to", tos)


toflip = c(1:7, 8, 10, 12:17)
cn[toflip] = paste("from",tos, "to", froms)[toflip]
fc_df[,-1][,toflip] = -fc_df[,-1][,toflip]

colnames(fc_df)[-1] = cn

###count files
setwd("~/zdrive/RNA_seq_data/Breast/")
count_files = c("All_5_BC_cell_Lines/v25/DESeq2_5BC_normalized_counts_v25.txt",
                "MCF10A_shRunx1/DESeq2/DESeq2_MCF10A_shRunx1_normalized_counts_v25.txt",
                "MCF10A-AT1_CA1_shRUNX2/DESeq2/MCF10A-AT1_shRUNX2/DESeq2_MCF10A-AT1_shRUNX2_normalized_counts_v25.txt",
                "MCF10A-AT1_CA1_shRUNX2/DESeq2/MCF10A-CA1_shRUNX2/2-way_comparison_EV_and_shRUNX2_only/DESeq2_MCF10A-CA1_shRUNX2_normalized_counts_v25.txt",
                "MCF7_shRunx1/DESeq2/DESeq2_MCF7_shRunx1_normalized_counts_v25.txt",
                "MDA-MB-231_Runx1/DESeq2/MDA_Runx1_3replicates_siRunx1_v21/DESeq2_3replicate_MDA_siRunx1_normalized_counts.txt")
names(count_files) = c("5BC", "MCF10A_shRunx1", "AT1_shRunx2", "CA1_shRunx2", "MCF7_shRunx1", "MDA231_siRunx1")
count_df = lapply(count_files, function(f)read.table(f, header = T))
#fix one with different header
count_df$MDA231_siRunx1= read.table(count_files["MDA231_siRunx1"], row.names = 1, header = T)

setwd(odir)
#convert id to stable
count_df = lapply(count_df, function(x){
    rownames(x) = sapply(strsplit(rownames(x), "\\."), function(str)str[[1]])
    x
})
#filter to match master
count_df = lapply(count_df, function(df){
    df = df[olaps_master$stable_gene_id,]
    rownames(df) = olaps_master$stable_gene_id
    df
})

### merge binding and DE
master = merge(olaps_master, fc_df, all.x = T)
master = merge(master, ref_dt[, c("stable_gene_id", "gene_type", "gene_name")])
o = c(ncol(master), ncol(master)-1, 1, 2:(ncol(master) -2 ))

master = master[,o, with = F]

library(openxlsx)
wb = createWorkbook()
options("openxlsx.numFmt" = "0.00")
s_head = createStyle(wrapText = T, valign = "top")

addWorksheet(wb, "binding and DE")
writeDataTable(wb, 1, master)
setColWidths(wb, 1, cols = 1:3, widths = 20)
setRowHeights(wb, 1, 1, 60)
addStyle(wb, 1, s_head, rows = 1, cols = 1:ncol(master))


for(i in seq_along(count_df)){
    s_name = paste(names(count_df)[i], "counts")
    addWorksheet(wb,s_name)
    writeDataTable(wb, s_name, count_df[[i]], rowNames = T)
    setColWidths(wb, s_name, cols = 1, widths = 20)
    setRowHeights(wb, s_name, 1, 60)
    addStyle(wb, s_name, s_head, rows = 1, cols = 1:(ncol(count_df[[i]])) + 1)
}

saveWorkbook(wb, file = "Runx1_binding_and_DE_for_selected_lncs_KMT2.xlsx", overwrite = T)
