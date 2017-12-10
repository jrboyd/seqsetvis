server_example_data = function(input, output, session, set_features_file, set_features_name, set_bigwig){
  #examples
  observeEvent(input$ExampleMCF7_bza, {
    showNotification(ui = "This data is for 4 histone marks from the MCF7 cell line treated with bezadoxifene.", duration = 10, id = "Note_ExMCF7_bza", type = "warning")
    #setting set_features_file, set_features_name, and set_bigwig is sufficient for valid setup
    set_features_file(paste0(bed_path, "/MCF7_bza_500ext.bed"))
    set_features_name("MCF7_bza_500ext")
    MCF7bza_bws = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs", pattern = "MCF7_bza_.+_FE.bw", full.names = T)
    names(MCF7bza_bws) = sub("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs/MCF7_bza_", "", MCF7bza_bws) %>%
      sub(pattern = "_FE.bw", replacement = "")
    example_bw = data.frame(filename = names(MCF7bza_bws), filepath = MCF7bza_bws, stringsAsFactors = F)
    set_bigwig(example_bw)
  })
  observeEvent(input$ExampleKasumi, {
    showNotification(ui = "This data is Dan Trombly's ChIPseq in Kasumi cell lines. For Kaleem.", duration = 10, id = "Note_ExKasumi", type = "warning")
    #setting set_features_file, set_features_name, and set_bigwig is sufficient for valid setup
    set_features_file(paste0(bed_path, "/Kasumi_bivalency.bed"))
    set_features_name("Kasumi_AE_bivalency")
    ex_bws = dir("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/aml/hg38/kasumi/DT_aml-eto_hg38/", pattern = "_FE.bw", full.names = T)
    names(ex_bws) = basename(ex_bws) %>% sub("Kasumi1_", "", .) %>%
      sub(pattern = "_pooled_FE.bw", replacement = "")
    example_bw = data.frame(filename = names(ex_bws), filepath = ex_bws, stringsAsFactors = F)
    set_bigwig(example_bw)
  })
  observeEvent(input$ExampleBivalency, {
    showNotification(ui = "Bivalency within 2kb of TSSes.", duration = 10, id = "Note_ExBiv", type = "warning")
    #setting set_features_file, set_features_name, and set_bigwig is sufficient for valid setup
    set_features_file(paste0(bed_path, "/TSS_serial_bivalency_2kb_ext.bed"))
    set_features_name("TSS_2kb_ext")
    ex_bws = c(
      "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K4ME3_pooled/H7_H3K4ME3_pooled_FE.bw",
      "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/H7_H3K27ME3_pooled/H7_H3K27ME3_pooled_FE.bw",
      "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/patients_GOOD1-H3K4ME3_pooled/patients_GOOD1-H3K4ME3_pooled_FE.bw",
      "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/patients_GOOD1-H3K27ME3_pooled/patients_GOOD1-H3K27ME3_pooled_FE.bw",
      "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/patients_POOR5-H3K4ME3_pooled/patients_POOR5-H3K4ME3_pooled_FE.bw",
      "/slipstream/galaxy/uploads/working/qc_framework/output_bivalency_redo_patients_H7/patients_POOR5-H3K27ME3_pooled/patients_POOR5-H3K27ME3_pooled_FE.bw",
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs/MCF7_ctrl_H3K4ME3_FE.bw",
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF7_drug_treatments_pooled_inputs/MCF7_ctrl_H3K27ME3_FE.bw",
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF10A_drug_treatments_pooled_inputs/MCF10A_ctrl_H3K4ME3_FE.bw",
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/breast/MCF10A_drug_treatments_pooled_inputs/MCF10A_ctrl_H3K27ME3_FE.bw"
    )
    names(ex_bws) = basename(ex_bws) %>% sub("Kasumi1_", "", .) %>%
      sub(pattern = "_pooled_FE.bw", replacement = "") %>%
      sub(pattern = "ctrl_", replacement = "") %>%
      sub(pattern = "OOR", replacement = "") %>%
      sub(pattern = "OOD", replacement = "") %>%
      sub(pattern = "_FE.bw", replacement = "") %>%
      sub(pattern = "patients_", replacement = "") %>%
      sub(pattern = "-", replacement = "_")
    example_bw = data.frame(filename = names(ex_bws), filepath = ex_bws, stringsAsFactors = F)
    set_bigwig(example_bw)
  })
  
  observeEvent(input$ExampleRunxAndCTCF, {
    showNotification(ui = "Bivalency within 2kb of TSSes.", duration = 10, id = "Note_ExBiv", type = "warning")
    #setting set_features_file, set_features_name, and set_bigwig is sufficient for valid setup
    set_features_file(paste0(bed_path, "/CTCF_across_10A_with_MK_MDA231_RUNX2.bed"))
    set_features_name("CTCF_and_RUNX")
    ex_bws = c(
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks///breast/10A_progression/AF-MCF10A_RUNX1_pooled_FE.bw",
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks///breast/10A_progression/AF-MCF10AT1_RUNX1_pooled_FE.bw",
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks///breast/10A_progression/AF-MCF10CA1_RUNX1_pooled_FE.bw",
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks///breast/10A_progression/MCF10A_CTCF_pooled_FE.bw",
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks///breast/10A_progression/MCF10AT1_CTCF_pooled_FE.bw",
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks///breast/10A_progression/MCF10CA1_CTCF_pooled_FE.bw",
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks///breast/MDA231_MK_runx/MDA231_Runx1_pooled_FE.bw",
      "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks///breast/MDA231_MK_runx/MDA231_Runx2_pooled_FE.bw"
    )
    names(ex_bws) = basename(ex_bws) %>% 
      gsub(pattern = "-", replacement = "_") %>%
      sub("Kasumi1_", "", .) %>%
      sub(pattern = "_pooled_FE.bw", replacement = "") %>%
      sub(pattern = "RUNX", replacement = "Runx") %>%
      sub(pattern = "AF_", replacement = "") %>%
      sub(pattern = "ctrl_", replacement = "") %>%
      sub(pattern = "MCF10CA1", replacement = "CA1") %>%
      sub(pattern = "MCF10AT1", replacement = "at1") %>%
      sub(pattern = "OOD", replacement = "") %>%
      sub(pattern = "_FE.bw", replacement = "") %>%
      sub(pattern = "patients_", replacement = "")
      
    example_bw = data.frame(filename = names(ex_bws), filepath = ex_bws, stringsAsFactors = F)
    set_bigwig(example_bw)
  })
  
  observeEvent(input$ExampleMCF10A, {
    showNotification(ui = "MCF10A histone marks (Terri). CTCF and RUNX1 (Andy).", duration = 10, id = "Note_MCF10A", type = "warning")
    #setting set_features_file, set_features_name, and set_bigwig is sufficient for valid setup
    set_features_file("/slipstream/home/ajfritz/ShinyData/NSvsC4_BETA_intersectRed.bed")
    set_features_name("MCF10A_CTCF_and_RUNX")
    ex_bws = c(
      "/slipstream/home/ajfritz/ShinyData/MCF10A_CTCF_pooled_FE.bw",     
      "/slipstream/home/ajfritz/ShinyData/MCF10A_H3K27ME3_pooled_FE.bw",  
      "/slipstream/home/ajfritz/ShinyData/MCF10A_H3K4ME3_pooled_FE.bw",
      "/slipstream/home/ajfritz/ShinyData/MCF10A_H4K20ME3_pooled_FE.bw",
      "/slipstream/home/ajfritz/ShinyData/MCF10A_H3K27AC_pooled_FE.bw",  
      "/slipstream/home/ajfritz/ShinyData/MCF10A_H3K4AC_pooled_FE.bw",    
      "/slipstream/home/ajfritz/ShinyData/MCF10A_H4K12AC_pooled_FE.bw",  
      "/slipstream/home/ajfritz/ShinyData/MCF10A_RUNX1_pooled_FE.bw"
    )
    names(ex_bws) = basename(ex_bws) %>% 
      gsub(pattern = "-", replacement = "_") %>%
      sub("Kasumi1_", "", .) %>%
      sub(pattern = "_pooled_FE.bw", replacement = "") %>%
      sub(pattern = "RUNX", replacement = "Runx") %>%
      sub(pattern = "AF_", replacement = "") %>%
      sub(pattern = "ctrl_", replacement = "") %>%
      sub(pattern = "MCF10CA1", replacement = "CA1") %>%
      sub(pattern = "MCF10AT1", replacement = "at1") %>%
      sub(pattern = "OOD", replacement = "") %>%
      sub(pattern = "_FE.bw", replacement = "") %>%
      sub(pattern = "patients_", replacement = "")
    
    example_bw = data.frame(filename = names(ex_bws), filepath = ex_bws, stringsAsFactors = F)
    set_bigwig(example_bw)
  })
  
  observeEvent(input$ExampleU937, {
    showNotification(ui = "U937 data. RUNX1 and IKZF1 overlap.", duration = 10, id = "Note_U937", type = "warning")
    #setting set_features_file, set_features_name, and set_bigwig is sufficient for valid setup
    set_features_file("/slipstream/home/joeboyd/ShinyData/IKZF1_Runx1_FE_10_olap.bed")
    set_features_name("U937_IKZF1_Runx1")
    ex_bws = c(
      "U937_IKZF1" =	"/slipstream/home//joeboyd/ShinyData//U937_data/U937_IK1h100_pooled_FE.bw",
      "U937_Runx1" =	"/slipstream/home//joeboyd/ShinyData//U937_data/U937_Runx1_pooled_FE.bw",
      "U937_H3K4me3" =	"/slipstream/home//joeboyd/ShinyData//U937_data/U937_H3K4me3_R1_FE.bw",
      "U937_H3K27me3" =	"/slipstream/home//joeboyd/ShinyData//U937_data/U937_H3K27me3_R1_FE.bw",
      "U937_H3K27ac" =	"/slipstream/home//joeboyd/ShinyData//U937_data/U937_H3K27ac_R1_FE.bw"
    )
    example_bw = data.frame(filename = names(ex_bws), filepath = ex_bws, stringsAsFactors = F)
    set_bigwig(example_bw)
  })
}