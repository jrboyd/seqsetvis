server_scatter_plot = function(input, output, session, 
                               get_features, 
                               get_profiles, 
                               set_selected){
  #the sampled and unsampled data being plotted in the scatterplot
  xy_plot_dt = reactiveVal({NULL})
  seed = reactiveVal(0)
  
  #the visible scatterplot data after sampling
  visible_dt = reactive({
    xy_val = xy_plot_dt()
    if(is.null(xy_val)) return(NULL)
    n_displayed = min(2000, nrow(xy_val))
    set.seed(seed())
    xy_val[sample(x = 1:nrow(xy_val), size = n_displayed)]
  })
  
  #the data selected by brushing
  selected_dt = observe({
    xy_dt = xy_plot_dt()
    if(is.null(xy_dt)) return(NULL)
    if(is.null(input$xy_brush)){
      sel_dt = xy_dt
    }else{
      sel_dt = brushedPoints(xy_dt, input$xy_brush)
    }
    set_selected(sel_dt)
  })
  
  #updates scatterplot when data to be plotted changes.
  output$xy_values = renderPlot({
    plotted_dt = visible_dt()
    if(is.null(plotted_dt)) return(NULL)
    #isolate because these trigger updates of xy_plot_dt() which already triggers reactivity
    x_variable = isolate(input$x_variable)
    y_variable = isolate(input$y_variable)
    ptype = input$RadioScatterplotType
    showHelp = input$CheckShowHelpers
    fixCoords = input$CheckFixCoords
    save(plotted_dt, showHelp, ptype, x_variable, y_variable, file = "last_plot.save")
    if(ptype == "standard"){
      plotted_dt[, max(xval, yval)]
      lim = c(0, plotted_dt[, max(xval, yval)])
      p = ggplot(plotted_dt) + 
        geom_point(aes(x = xval, y = yval, col = plotting_group)) +
        labs(x = x_variable, y = y_variable, title = "Max FE in regions") +
        ylim(lim) + xlim(lim)
      if(fixCoords) p = p + coord_fixed()
      if(showHelp){
        pos1 = .2 * max(lim)
        pos2 = .8 * max(lim)
        p = p + annotate("segment", x = 0, xend = 0, y = pos1, yend = pos2, arrow = arrow())
        p = p + annotate("label", x = 0, y = mean(lim), label = gsub(" ", "\n", paste(y_variable, "binding")), hjust = 0)
        p = p + annotate("segment", x = pos1, xend = pos2, y = 0, yend = 0, arrow = arrow())
        p = p + annotate("label", x = mean(lim), y = 0, label = paste(x_variable, "binding"), vjust = 0)
        p = p + annotate("label", x = 0, y = 0, label = "no\nbinding", hjust = 0, vjust = 0)
        p = p + annotate("segment", x = pos1, xend = pos2, y = pos1, yend = pos2, arrow = arrow())
        p = p + annotate("label", x = mean(lim), y = mean(lim), label = gsub(" ", "\n", paste("both binding")))
        p
      }
    }else if(ptype == "volcano"){
      # plotted_dt[, xvolcano := log2(max(yval, 1) / max(xval, 1)), by = id]
      # plotted_dt[, yvolcano := max(min(yval, xval), 1), by = id]
      xmax = plotted_dt[, max(abs(c(xvolcano)))]
      lim = c(-xmax, xmax)
      p = ggplot(plotted_dt) + 
        geom_point(aes(x = xvolcano, y = yvolcano, col = plotting_group)) +
        labs(x = paste("log2 fold-change of", y_variable, "over", x_variable), 
             y = paste("log2 min of", y_variable, "and", x_variable), title = "Max FE in regions") +
        xlim(lim)
      if(fixCoords) p = p + coord_fixed()
      if(showHelp){
        pos1 = .2 * max(lim)
        pos2 = .8 * max(lim)
        p = p + annotate("segment", y = 1, yend = 1, x = pos1, xend = pos2, arrow = arrow())
        p = p + annotate("label", y = 1, x = max(lim)/2, label = gsub(" ", "\n", paste(y_variable, "binding")), vjust = 0)
        p = p + annotate("segment", y = 1, yend = 1, x = -pos1, xend = -pos2, arrow = arrow())
        p = p + annotate("label", y = 1, x = -max(lim)/2, label = gsub(" ", "\n", paste(x_variable, "binding")), vjust = 0)
        p = p + annotate("label", x = 0, y = 1, label = "no\nbinding", hjust = .5, vjust = 0)
        ylim = range(plotted_dt$yvolcano)
        ypos1 = 1 + (max(ylim) - min(ylim)) * .2
        ypos2 = 1 + (max(ylim) - min(ylim)) * .8
        p = p + annotate("segment", x = 0, xend = 0, y = ypos1, yend = ypos2, arrow = arrow())
        p = p + annotate("label", x = mean(lim), y = mean(ylim), label = gsub(" ", "\n", paste("both binding")))
        p
      }
    }
    p
  })
  
  
  
  #set xy data when appropriate
  observe({
    if(is.null(input$GroupingType)) return(NULL)
    prof_dt = get_profiles()
    #TODO selector for this
    # samples_loaded = unique(prof_dt$sample)
    x_variable = input$x_variable
    y_variable = input$y_variable
    if(is.null(x_variable) || is.null(y_variable)) return(NULL)
    x = prof_dt[sample == x_variable]
    y = prof_dt[sample == y_variable]
    x_summit_val = x[, .(xval = FE[which(FE == max(FE, na.rm = T))[1]]), by = .(hit)]#[order(as.numeric(hit))]
    y_summit_val = y[, .(yval = FE[which(FE == max(FE, na.rm = T))[1]]), by = .(hit)]#[order(as.numeric(hit))]
    xy_val = merge(x_summit_val, y_summit_val, by = "hit")

    # xy_val$hit = as.integer(xy_val$hit)
    feat_gr = isolate(get_features())
    feat_dt = as.data.table(feat_gr)
    # feat_dt$id = as.integer(feat_dt$id)
    setkey(feat_dt, id)
    xy_val = cbind(xy_val, feat_dt[.(xy_val$hit),])
    
    if(input$GroupingType == "logical derived"){
      logical_groups = input$SelectLogicalGrouping
      
      switch(paste0("l", length(logical_groups)),
             l0 = {
               xy_val$plotting_group = "no plotting_group set"
             },
             l1 = {
               print(1)
               xy_val$plotting_group = "negative"
               xy_val[get(logical_groups[1]), plotting_group := logical_groups[1]]
             },
             l2 = {
               print(2)
               xy_val$plotting_group = "neither"
               xy_val[get(logical_groups[1]) & get(logical_groups[2]), plotting_group := "both"]
               xy_val[get(logical_groups[1]) & !get(logical_groups[2]), plotting_group := logical_groups[1]]
               xy_val[!get(logical_groups[1]) & get(logical_groups[2]), plotting_group := logical_groups[2]]
             },
             l3 = {
               print(3)
               xy_val$plotting_group = "none"
               xy_val[get(logical_groups[1]) & get(logical_groups[2]) & get(logical_groups[3]), plotting_group := "all"]
               xy_val[!get(logical_groups[1]) & get(logical_groups[2]) & get(logical_groups[3]), plotting_group := paste(logical_groups[2:3], collapse = " & ")]
               xy_val[get(logical_groups[1]) & !get(logical_groups[2]) & get(logical_groups[3]), plotting_group := paste(logical_groups[c(1,3)], collapse = " & ")]
               xy_val[get(logical_groups[1]) & get(logical_groups[2]) & !get(logical_groups[3]), plotting_group := paste(logical_groups[1:2], collapse = " & ")]
               xy_val[!get(logical_groups[1]) & !get(logical_groups[2]) & get(logical_groups[3]), plotting_group := logical_groups[3]]
               xy_val[get(logical_groups[1]) & !get(logical_groups[2]) & !get(logical_groups[3]), plotting_group := logical_groups[1]]
               xy_val[!get(logical_groups[1]) & get(logical_groups[2]) & !get(logical_groups[3]), plotting_group := logical_groups[2]]
             },
             {
               print("too many groups")
               xy_val$plotting_group = "too many groups"
             }
      )
    }else if(input$GroupingType == "predefined"){
      fgrp_name = input$SelectFactorGrouping
      if(fgrp_name == "chromosomes"){
        fgrp = xy_val$seqnames
      }else{
        fgrp = xy_val[[fgrp_name]]
      }
      xy_val$plotting_group = fgrp
      
    }else{
      stop(paste("bad GroupingType", input$GroupingType))
    }
    xy_val[, xvolcano := log2(max(yval, 1) / max(xval, 1)), by = id]
    xy_val[, yvolcano := log2(max(min(yval, xval), 1)), by = id]
    xy_plot_dt(xy_val)
  })
  
  #generates UI appropriate to metadata columns in get_features
  output$GroupingsAvailable = renderUI({
    fgr = get_features()
    chrms = unique(seqnames(fgr))
    mdat = elementMetadata(fgr)
    if(!is.null(mdat$id)) mdat$id = NULL
    if(ncol(mdat) > 0){
      col_classes = sapply(1:ncol(mdat), function(i)class(mdat[,i]))
      #groups composed of just T and F can be compared to create new sets
      logical_groups = sort(colnames(mdat)[which(col_classes == "logical")])
      #groups of factors can only be analyzed individually
      factor_names = colnames(mdat)[which(col_classes != "logical")]
      names(factor_names) = factor_names
      factor_groups = list(chromosomes = chrms)
      factor_groups = append(factor_groups, lapply(factor_names, function(x){
        unique(mdat[[x]])
      }))
    }else{
      logical_groups = character()
      factor_groups = list(chromosomes = chrms)
    }
    factor_groups = lapply(factor_groups, as.character)
    if(length(logical_groups) > 0){
      #create a conditional selector
      #the factor grouping chromosome (GRanges seqnames) is always assumed present
      #no radio button will be shown if logicals aren't present
      tagList(
        radioButtons(inputId = "GroupingType", label = "Groupings Available", choices = c("logical derived", "predefined"), selected = "logical derived"),
        conditionalPanel(
          condition = "input.GroupingType == 'logical derived'",
          selectInput("SelectLogicalGrouping", label = "Select Logical Groups", choices = logical_groups, multiple = T, selectize = T, selected = logical_groups[1:max(1, min(2, length(logical_groups)))])
        ),
        conditionalPanel(
          condition = "input.GroupingType == 'predefined'",
          selectInput("SelectFactorGrouping", label = "Select Factor Group", choices = names(factor_groups), selected = names(factor_groups)[1], multiple = F, selectize = F)
        )
      )
    }else{
      tagList(
        (radioButtons(inputId = "GroupingType", label = "Groupings Available", choices = c("predefined"), selected = "predefined")),
        selectInput("SelectFactorGrouping", label = "Select Factor Group", 
                    choices = names(factor_groups), selected = names(factor_groups)[1], 
                    multiple = F, selectize = F)
      )
      
    }
    
  })
}
