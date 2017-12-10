#all of the non-reactive server elements
library(shiny)
library(shinyjs)
library(shinyFiles)
library(GenomicRanges)
library(ggplot2)
library(data.table)
library(magrittr)
source("functions_intersect.R")
source("source_gg_venneuler.R")
source("jrb_gg_vennDiagram.R")

bed_path = "~/ShinyApps/shiny_peak_data/beds"
dir.create(bed_path, showWarnings = F)
options(shiny.maxRequestSize=50*1024^2)

# Return the UI for a modal dialog with data selection input. If 'failed' is
# TRUE, then display a message that the previous value was invalid.
dataModal <- function(sets, failed = FALSE) {
  sets_html = HTML(paste(sapply(sets$selected[1], function(x){
    as.character(textInput("TxtRename", label = paste("rename", x), value = x))
  }), collapse = "\n"))
  modalDialog(
    sets_html,
    span('(Please rename the selected sample)'),
    if (failed)
      div(tags$b("One or more names no longer unique!", style = "color: red;")),
    
    footer = tagList(
      modalButton("Cancel"),
      actionButton("BtnConfirmRename", "Confirm")
    )
  )
}



shinyFiles2load = function(shinyF, roots){
  root_path = roots[shinyF$root]
  rel_path = paste0(unlist(shinyF$files), collapse = "/")
  file_path = paste0(root_path, "/", rel_path)
  return(file_path)
}

shinyFiles2save = function(shinyF, roots){
  root_path = roots[shinyF$root]
  rel_path = paste0(unlist(shinyF$name), collapse = "/")
  file_path = paste0(root_path, "/", rel_path)
  return(file_path)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#assumed colnames for MACS2 peak files
peak_cn = c("seqnames", "start", "end", "id", "score", "strand", "FE", "p-value", "q-value", "summit_pos")

get_col_classes.df = function(df){
  sapply(1:ncol(df), function(i)class(df[[i]]))
}
get_col_classes = function(file, skipFirst = F){
  df = read.table(file, stringsAsFactors = F, header = F, nrows = 10, skip = ifelse(skipFirst, 1, 0))
  get_col_classes.df(df)
}
test_for_header = function(peak_file){
  withFirst = get_col_classes(peak_file, skipFirst = F)
  skipFirst = get_col_classes(peak_file, skipFirst = T)
  #definitely no header
  if(all(withFirst == skipFirst)){
    return(F)
  }
  if(all(withFirst == "character")){
    return(T)
  }
  warning(paste("can't determine if", peak_file, "has a header, guess not."))
  return(F)
}
load_peak_wValidation = function(peak_file, with_notes = F){
  has_header = test_for_header(peak_file)
  df = read.table(peak_file, stringsAsFactors = F, header = has_header)
  col_classes = get_col_classes.df(df)
  if(!has_header){
    if(ncol(df) == length(peak_cn)){
      if(with_notes){
        showNotification("assuming file is narrowPeak.", type = "warning")
      }else{
        print("assuming file is narrowPeak.")
      }
      colnames(df) = peak_cn  
    }else{
      if(with_notes){
        showNotification("file not narrowPeak, loading as minimal bed file.", type = "warning")
      }else{
        print("file not narrowPeak, loading as minimal bed file.")
      }
      bed_cn = peak_cn[1:4]
      nc = min(ncol(df), length(bed_cn))
      colnames(df)[1:nc] = peak_cn[1:nc]
    }
  }else{#try to make colnames GRanges compatible
    if(all(col_classes[1:5] == c("character", "integer", "integer", "integer", "character"))){
      if(with_notes){
        showNotification("file looks like saved GRanges.", type = "warning")
      }else{
        print("file looks like saved GRanges.")
      }
      colnames(df)[1:5] = c("seqnames", "start", "end", "width", "strand")
    }else if(all(col_classes[1:ncol(df)] == c("character", "integer", "integer", 
                                              "character", "integer", "character", 
                                              "numeric", "numeric", "numeric", "integer")[1:ncol(df)])){
      if(with_notes){
        showNotification("file looks like bed or encode peak", type = "warning")
      }else{
        print("file looks like bed or encode peak")
      }
    }else{
      colnames(df)[1:3] = c("seqnames", "start", "end")
      if(with_notes){
        showNotification("forced to assume first 3 columns are minimal bed, might break.", type = "warning")
      }else{
        print("forced to assume first 3 columns are minimal bed, might break.")
      }
    }
  }
  print(head(df))
  print(paste0(nrow(df), " total rows..."))
  invisible(df)
}

#setup root paths
user_roots = dir("/slipstream/home/", full.names = T) %>% paste0(. , "/ShinyData")
user_roots = subset(user_roots, dir.exists(user_roots))
names(user_roots) = dirname(user_roots) %>% basename()
qcframework_load <<- dir("/slipstream/galaxy/uploads/working/qc_framework", pattern = "^output", full.names = T)
names(qcframework_load) <- basename(qcframework_load)
roots_load = c(user_roots, qcframework_load)
roots_output =  c("intersectR" = bed_path, user_roots)

create_metaDF_empty = function(){
  df = data.frame(id_name = character(), 
                  data_frame = I(list()), 
                  file_path = character(), 
                  display_name = character(), stringsAsFactors = F)
  return(df)
}
create_metaDF_row = function(id_name, data_frame, file_path, display_name){
  if(class(data_frame) != "list") data_frame = list(data_frame)
  df = data.frame(id_name = id_name, 
                  data_frame = I(data_frame), 
                  file_path = file_path, 
                  display_name = display_name, 
                  stringsAsFactors = F,
                  row.names = id_name)
  return(df)
}