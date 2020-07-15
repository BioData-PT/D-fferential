source("mod_filter_table.R")

########################################
# edgeR options UI
mod_edgeRClassicUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    h4("edgeR Classic Options"),
    fluidRow(
      column(4, textInput(ns("txt_conditionA"), "Condition A", value = "ConditionA")),
      column(8, selectizeInput(ns("sel_conditionA"), "Columns", choices=NULL, multiple = TRUE))
    ),
    fluidRow(
      column(4, textInput(ns("txt_conditionB"), "Condition B", value = "ConditionB")),
      column(8, selectizeInput(ns("sel_conditionB"), "Columns", choices=NULL, multiple = TRUE))
    ),
    fluidRow(
      column(6, numericInput(ns("num_fdr"), "FDR", value = 0.05)),
      column(6, numericInput(ns("num_logfc"), "LogFC", value = 1))
    )
  )
}

########################################
# edgeR module
# TODO: we are getting all inputs from the tab_de module. Idealy we would only need a few...
mod_edgeRClassicServer <- function(input, output, session, dataframe, de.input) {
  
  # fill in the column choices when dataframe changes
  observe({
    df <- dataframe()
    cols <- colnames(df)
    
    # preserve previously selected columns
    isolate({
      selA <- input$sel_conditionA
      selB <- input$sel_conditionB
      
      updateSelectizeInput(session, "sel_conditionA", choices=cols, selected = intersect(cols, selA))
      updateSelectizeInput(session, "sel_conditionB", choices=cols, selected = intersect(cols, selB))
    })
  })
  
  # remove option from B when added to A
  observe({
    colsA <- input$sel_conditionA
    
    isolate({
      df <- dataframe()
      cols <- colnames(df)
      cols <- setdiff(cols, colsA)
      
      selected <- input$sel_conditionB
      
      updateSelectizeInput(session, "sel_conditionB", choices=cols, selected = selected)
    })
  })
  
  # remove option from A when added to B
  observe({
    colsB <- input$sel_conditionB
    
    isolate({
      df <- dataframe()
      cols <- colnames(df)
      cols <- setdiff(cols, colsB)
      
      selected <- input$sel_conditionA
      
      updateSelectizeInput(session, "sel_conditionA", choices=cols, selected = selected)
    })
  })

  # build the expression matrix
  expr.matrix <- reactive({
    colsA <- input$sel_conditionA
    colsB <- input$sel_conditionB
    df <- dataframe()
    
    validate(need(length(colsA) > 0, message = FALSE),
             need(length(colsB) > 0, message = FALSE))
    
    return (cbind(df[, match(colsA, colnames(df)), drop=FALSE ],
                  df[, match(colsB, colnames(df)), drop=FALSE ]))
  })
  
  condition.factor <- reactive({
    colsA <- input$sel_conditionA
    colsB <- input$sel_conditionB
    nameA <- input$txt_conditionA
    nameB <- input$txt_conditionB
    
    validate(need(length(colsA) > 1, message = FALSE),
             need(length(colsB) > 1, message = FALSE))
    
    cfact <- factor(c(rep(nameA, length(colsA)), rep(nameB, length(colsB))), 
                    levels = c(nameA, nameB))
    
    return(cfact)
  })
  
  # do edgeR when button is clicked
  result <- reactive({
    counts <- expr.matrix()
    conditions <- condition.factor()
    
    # do edgeR
    withProgress(message = "Running edgeR...", value = 0, {
      n = 5
      
      y <- DGEList(counts=counts, group = conditions)
      
      incProgress(1/n, detail = "Calculating norm factors...")
      y <- calcNormFactors(y)
      
      #--------------------------------------------------------------------
      # aggs: perform quality-control assessment 
      
      incProgress(1/n, detail = "Plot MDS...")
      mds <- plotMDS( x = y, plot = FALSE) # MDS data 
    
      colData <- data.frame( "condition" = conditions) # create meta tbl df
      rownames(colData) <- colnames(counts) # add row names to df
      
      # create a data frame with the MDS (PCoA) coordinates
      mds_coord <- data.frame("Group" = names(mds$x), 
                               #"Samples" = rownames(colData),
                               "Condition" = colData$condition,
                               "x_axis" = mds$x, 
                               "y_axis" = mds$y)
      
      # use ggplot to plot MDS (PCoA) coordinates - condition var
      dim_plot <- ggplot(mds_coord,
                         aes(x = x_axis, 
                             y = y_axis, 
                             label = Group)) + 
        geom_point(aes(color = Condition), 
                   size = 3) + 
        ggrepel::geom_text_repel(size = 3.5) +
        xlab("Leading logFC dim 1") + 
        ylab("Leading logFC dim 2")
      
      incProgress(1/n, detail = "Calculating CPM to heatmap...")
      logcpm <- cpm( y = y, log = TRUE )
      norm_cpm <- cpm( y = y )
      
      corr_mtx <- cor( logcpm ) # mtx and correlation mtx of logCPM trans
  
      incProgress(1/n, detail = "Plot sample-to-sample heatmap...")
      heat_plot <- heatmaply::ggheatmap( corr_mtx, row_dend_left = TRUE, 
                                         col_side_colors = colData$condition ) # sample-to-sample heatmap
      
      #--------------------------------------------------------------------
      
      incProgress(1/n, detail = "Estimating dispersions...")
      y <- estimateDisp(y)
      
      incProgress(1/n, detail = "Performing exact test...")
      et <- exactTest(y)
      
      incProgress(1/n, detail = "Extracting results...")
      topgenes <- topTags(et, n=dim(counts)[[1]])
    })
    
    # save result
    # n <- length(names(results)) + 1
    
    #if (nchar(input$text_name) > 0) {
    id <- de.input$txt_name
    #} else {
    #  id <- paste("Differential expression", n)
    #}
    
    result <- list(
      id = id,
      error = NA,
      method = "edgeR classic",
      samples = colnames(counts),
      conditions = conditions,
      fdr = input$num_fdr,
      lfc = input$num_logfc,
      table = topgenes$table,
      plot_dim = dim_plot, # aggs
      plot_heat = heat_plot, # aggs 
      corr_mtx = corr_mtx, # aggs
      ge_mtx_trans = logcpm, # aggs
      ge_mtx_norm = norm_cpm, # aggs
      genes.up = rownames(topgenes$table)[ which(topgenes$table$FDR < input$num_fdr & topgenes$table$logFC >= input$num_logfc) ],
      genes.down = rownames(topgenes$table)[ which(topgenes$table$FDR < input$num_fdr & topgenes$table$logFC <= -input$num_logfc) ],
      genes.diff = rownames(topgenes$table)[ which(topgenes$table$FDR < input$num_fdr & abs(topgenes$table$logFC) >= input$num_logfc) ]
    )
    
    return(result)
  })
  
  ready <- reactive({
    colsA <- input$sel_conditionA
    colsB <- input$sel_conditionB
    
    if (length(colsA) < 2 || length(colsB) < 2)
      return(FALSE)
    else {
      return(TRUE)
    }
  })
  
  return(list(
    result = result,
    ready = ready
  ))
}


########################################
# edgeR options UI
mod_deseq2UI <- function(id) {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    h4("DESeq2 Options"),
    fluidRow(
      column(4, textInput(ns("txt_conditionA"), "Condition A", value = "ConditionA")),
      column(8, selectizeInput(ns("sel_conditionA"), "Columns", choices=NULL, multiple = TRUE))
    ),
    fluidRow(
      column(4, textInput(ns("txt_conditionB"), "Condition B", value = "ConditionB")),
      column(8, selectizeInput(ns("sel_conditionB"), "Columns", choices=NULL, multiple = TRUE))
    ),
    fluidRow(
      column(6, numericInput(ns("num_fdr"), "FDR", value = 0.05)),
      column(6, numericInput(ns("num_logfc"), "LogFC", value = 1))
    )
  )
}

########################################
# deseq2 module
# TODO: we are getting all inputs from the tab_de module. Idealy we would only need a few...
mod_deseq2Server <- function(input, output, session, dataframe, de.input) {
  
  # fill in the column choices when dataframe changes
  observe({
    df <- dataframe()
    cols <- colnames(df)
    
    # preserve previously selected columns
    isolate({
      selA <- input$sel_conditionA
      selB <- input$sel_conditionB
      
      updateSelectizeInput(session, "sel_conditionA", choices=cols, selected = intersect(cols, selA))
      updateSelectizeInput(session, "sel_conditionB", choices=cols, selected = intersect(cols, selB))
    })
  })
  
  # remove option from B when added to A
  observe({
    colsA <- input$sel_conditionA
    
    isolate({
      df <- dataframe()
      cols <- colnames(df)
      cols <- setdiff(cols, colsA)
      
      selected <- input$sel_conditionB
      
      updateSelectizeInput(session, "sel_conditionB", choices=cols, selected = selected)
    })
  })
  
  # remove option from A when added to B
  observe({
    colsB <- input$sel_conditionB
    
    isolate({
      df <- dataframe()
      cols <- colnames(df)
      cols <- setdiff(cols, colsB)
      
      selected <- input$sel_conditionA
      
      updateSelectizeInput(session, "sel_conditionA", choices=cols, selected = selected)
    })
  })
  
  # build the expression matrix
  expr.matrix <- reactive({
    colsA <- input$sel_conditionA
    colsB <- input$sel_conditionB
    df <- dataframe()
    
    validate(need(length(colsA) > 0, message = FALSE),
             need(length(colsB) > 0, message = FALSE))
    
    return (cbind(df[, match(colsA, colnames(df)), drop=FALSE ],
                  df[, match(colsB, colnames(df)), drop=FALSE ]))
  })
  
  condition.factor <- reactive({
    colsA <- input$sel_conditionA
    colsB <- input$sel_conditionB
    nameA <- input$txt_conditionA
    nameB <- input$txt_conditionB
    
    validate(need(length(colsA) > 1, message = FALSE),
             need(length(colsB) > 1, message = FALSE))
    
    cfact <- factor(c(rep(nameA, length(colsA)), rep(nameB, length(colsB))), 
                    levels = c(nameA, nameB))
    
    return(cfact)
  })
  
  
  # do edgeR when button is clicked
  result <- reactive({
    counts <- expr.matrix()
    conditions <- condition.factor()
    
    # do edgeR
    withProgress(message = "Running DESeq2...", {
      colData <- data.frame(colnames(counts),
                            condition=conditions,
                            row.names=1)
      
      dds <- DESeqDataSetFromMatrix(countData = counts, 
                                    colData = colData, 
                                    design = ~condition)
      
      #--------------------------------------------------------------------
      # aggs: perform quality-control assessment 

      r_dds <- rlog( object = dds, blind = TRUE ) # rlog transformation - it's blind to the exp design
      pcaCoordData <- plotPCA( object = r_dds, 
                      intgroup = "condition", 
                      returnData = TRUE ) # pca data by condition to plot below
      varPCs <- round(100 * attr(pcaCoordData, "percentVar") ) # get PCs variance
      
      dim_plot <- ggplot(pcaCoordData,
                         aes(x = PC1, 
                             y = PC2, 
                             label = name)) + 
        geom_point(aes(color = group), 
                   size = 3) + 
        ggrepel::geom_text_repel(size = 3.5) +
        xlab(paste0("PC1 (",varPCs[1],"% explained variance)")) +
        ylab(paste0("PC2 (",varPCs[2],"% explained variance)")) 
      

      corr_mtx <- cor( assay( r_dds ) ) # mtx and correlation mtx of rlog trans
      heat_plot <- heatmaply::ggheatmap( corr_mtx, #row_dend_left = TRUE, 
                              col_side_colors = colData$condition ) # sample-to-sample heatmap
      
      
      ## normalize counts and saved it 
      dds_size_fct <- estimateSizeFactors( dds )
      norm_counts <- counts( dds_size_fct, normalized = TRUE )
      
      #--------------------------------------------------------------------
      
      dds <- DESeq(dds)
      res <- results(dds)
    })
    
    # save result
    # n <- length(names(results)) + 1
    
    #if (nchar(input$text_name) > 0) {
    id <- de.input$txt_name
    #} else {
    #  id <- paste("Differential expression", n)
    #}
    
    tab <- as.data.frame(res)
    
    result <- list(
      id = id,
      error = NA,
      method = "DESeq2",
      samples = colnames(counts),
      conditions = conditions,
      fdr = input$num_fdr,
      lfc = input$num_logfc,
      table = tab,
      plot_dim = dim_plot, # aggs
      plot_heat = heat_plot, # aggs
      corr_mtx = corr_mtx, # aggs
      ge_mtx_trans = assay(r_dds), # aggs
      ge_mtx_norm = norm_counts, # aggs
      genes.up = rownames(tab)[ which(tab$padj < input$num_fdr & tab$log2FoldChange >= input$num_logfc) ],
      genes.down = rownames(tab)[ which(tab$padj < input$num_fdr & tab$log2FoldChange <= -input$num_logfc) ],
      genes.diff = rownames(tab)[ which(tab$padj < input$num_fdr & abs(tab$log2FoldChange) >= input$num_logfc) ]
    )
    
    #results[[ id ]] <- result
    return(result)
  })
  
  ready <- reactive({
    colsA <- input$sel_conditionA
    colsB <- input$sel_conditionB
    
    if (length(colsA) < 2 || length(colsB) < 2)
      return(FALSE)
    else {
      return(TRUE)
    }
  })
  
  return(list(
    result = result,
    ready = ready
  ))
}


#' UI function for the de tab
tab_deUI <- function(id) {
  ns <- NS(id)
  
  de.options <- tagList(
    h4("Differential Expression Options"),
    selectInput(ns("sel_method"), "Method", 
                choices=c("edgeR classic" = "edger_classic", 
                          #"edgeR GLM" = "edgeR_glm", # aggs: uncomment it; it's not implemented yet. 
                          "DESeq2" = "deseq2"), selected = "edgeR classic"),
    hr(),
    conditionalPanel(condition = paste0("input['", ns("sel_method"), "']", " == 'edger_classic'"), 
                     mod_edgeRClassicUI(ns("mod_edger_classic"))),
    conditionalPanel(condition = paste0("input['", ns("sel_method"), "']", " == 'deseq2'"), 
                     mod_deseq2UI(ns("mod_deseq2"))),
    hr(),
    textInput(ns("txt_name"), "Name for this analysis"),
    actionButton(ns("btnGo"), "Go!", width="100%")
  )

  edger.result.ui <- tagList(
    htmlOutput(ns("html_res_params")),
    hr(),
    tabsetPanel(
      #---------------------------------------------------------------------------------------------------------------------
      # aggs: add QC plots tab
      tabPanel("QC Plots", 
               plotOutput(ns("dim_plot")),
               helpText("MDS (in the case of EdgeR) or PCA (in the case of DESeq2) plot", 
                        "showing the variability among groups. In the case of MDS", 
                        "the plot represents the leading log-fold-change between", 
                        "pair of samples based on the top 500 genes normalized.", 
                        "In the case of PCA, the plot represents the top 500 most",
                        "variable genes rlog transformed."), 
               hr(),
               plotOutput(ns("heat_plot")),
               helpText("Heatmap plot highlighting the sample-to-sample variability", 
                        "based on the Pearson correlation. In the case of edgeR,",
                        "the logCPM transformed counts were used to determine the correlations", 
                        "between sample-to-sample. In the case of DESeq2, the rlog transformed", 
                        "counts were used to determine the Pearson correlations between sample-to-sample.")
      ),
      #---------------------------------------------------------------------------------------------------------------------
      tabPanel("DE Plots", 
               # inputPanel(
               #   radioButtons(ns("radio_volcano_type"), label = "Plot type", )
               # ),
               plotOutput(ns("volcano_plot")),
               helpText("Volcano plot, showing the relationship between log2 fold-change",
                        "and evidence of differential expression."),
               hr(),
               plotOutput(ns("plot_ma")),
               helpText("MA plot, showing the relationship between mean expression",
                        "and log2 fold-change.")),
      #---------------------------------------------------------------------------------------------------------------------
      # aggs: add Downstream plots tab
      tabPanel("Downstream Plots", 
               plotOutput(ns("heat_de_plot")),
               helpText("Heatmap highlighting the clustering of genes (rows), by regulation type,", 
                        "and samples (columns), by condition, of genes differentially expressed.", 
                        "Only the top 25 genes up- and downregulated are shown (if available,",
                        "otherwise <25 or none, depending on the no. of DGE obtained).",
                        "In the case of edgeR, it uses the log2 CPM counts normalized and scaled.", 
                        "In the case of DESeq2, it uses the counts rlog transformed and scaled."), 
               hr(),
               plotOutput(ns("dot_dge_up_plot")),
               helpText("Dot plot highlighting the top upregulated genes that were found among the", 
                        "genes differentially expressed. The top upregulated genes means the genes", 
                        "with the highest log fold-change values.", 
                        "In the case of edgeR, the counts highlighted in the y-axis are normalized in CPM.", 
                        "In the case of DESeq2, it uses the DESeq2 normalization method.", 
                        "The y-axis appears log10 transformed.",
                        "In the interactive tab, the genes with 0 counts here were substituted by 1e-06",  
                        "in order to plot the interactive dot plot."), 
               hr(),
               plotOutput(ns("dot_dge_down_plot")),
               helpText("Dot plot highlighting the top downregulated genes that were found among the", 
                        "genes differentially expressed. The top downregulated genes means the genes", 
                        "with the lowest log fold-change values.", 
                        "In the case of edgeR, the counts highlighted in the y-axis are normalized in CPM.", 
                        "In the case of DESeq2, it uses the DESeq2 normalization method.", 
                        "The y-axis appears log10 transformed.",
                        "In the interactive tab, the genes with 0 counts here were substituted by 1e-06",  
                        "in order to plot the interactive dot plot.")),
      
      #---------------------------------------------------------------------------------------------------------------------
      tabPanel("Interactive Plots", 
               # inputPanel(
               #   radioButtons(ns("radio_volcano_type"), label = "Plot type", )
               # ),
               #ggvisOutput(ns("ggvis_volcano")),
               inputPanel(radioButtons(ns("radio_plotly_type"), "Plot type", 
                                       choices=c("MDS/PCA", "Heatmap - QC", 
                                                 "Volcano plot", "MA-plot", 
                                                 "Heatmap - DE",
                                                 "Dotplot - up", "Dotplot - down"))), # aggs: add QC plots
               plotlyOutput(ns("plotly_plot"))
               # hr(),
               # ggvisOutput(ns("ggvis_ma")),
               # helpText("MA plot, showing the relationship between mean expression",
               #          "and log2 fold-change.")),
      ),
      tabPanel("Table", 
               dataTableOutput(ns("result_table")),
               hr(),
               downloadButton(ns("download_table"), "Download (.csv)"))
      )
    )
  
  main.panel.ui <- tagList(
    fluidRow(
      column(8, selectInput(ns("sel_result"), label = "Show Result", choices = NULL, width = "100%")),
      column(4, actionButton(ns("button_delete"), label = "Delete this analysis"))),
    edger.result.ui,
    verbatimTextOutput(ns("debug"))
  )
  
  # sidebarLayout(
  #   sidebarPanel(de.options),
  #   mainPanel(main.panel.ui)
  # )
  fluidRow(
    box(de.options, width=4),
    conditionalPanel(condition = "true",
                     box(main.panel.ui, width=8))
  )
}

#' Server function for statistics module
#' 
#' @param cmatrix Counts matrix.
#' 
#' @return A dataframe as a reactive value.
tab_deServer <- function(input, output, session, sessionData) {

  dataframe <- sessionData$dataframe
  
  edgeR_classic <- callModule(mod_edgeRClassicServer, id = "mod_edger_classic", dataframe, input)
  deseq2 <- callModule(mod_deseq2Server, id = "mod_deseq2", dataframe, input)
  
  # start with the Go button disabled
  shinyjs::disable("btnGo")
  
  # this is where we store the analysis results
  results <- reactiveValues(data=list())
  
  # reset results when data changes
  observe({
    dataframe()
    
    print("dataframe changed")
    
    # results$data <- list()
  })
  
  # enable the Go button when required info is entered
  observe({
    ready <- switch(input$sel_method,
                    "edger_classic"=edgeR_classic$ready(),
                    "edgeR_glm"=FALSE,
                    "deseq2"=deseq2$ready())
    
    if (ready == FALSE)
      shinyjs::disable("btnGo")
    else {
      shinyjs::enable("btnGo")
    }
  })

  # observe new results and update the select input and name for next analysis  
  observe({
    print("updating results")
    
    updateTextInput(session, "txt_name", value = paste("Differential expression", length(names(results$data)) + 1))
    updateSelectInput(session, "sel_result", choices = names(results$data), selected = rev(names(results$data))[1])
  })
  
  # get the selected result
  result <- reactive({
    req(input$sel_result)
    
    validate(need(length(results$data) > 0, message = FALSE))
    
    results$data[[ input$sel_result ]]
  })
  
  # 
  result.table <- reactive({
    tab <- result()$tab
    
    tab$Result <- "No Change"
    tab[ result()$genes.up, ]$Result <- "Up"
    tab[ result()$genes.down, ]$Result <- "Down"

    tab
  })

  output$download_table <- downloadHandler(
    filename = function() {
      paste0("d-fferential-", Sys.Date(), "-", result()$id, ".csv")
    },
    content = function(file) {
      write.csv(result.table(), file)
    }
  )
  
  # show the results table
  output$result_table <- renderDataTable({
    tab <- result.table()
    
    DT::datatable(tab, options = list(scrollX = TRUE), filter="top")
  })
  
  # summary of the analysis
  output$html_res_params <- renderText({
    res <- result()
    
    paste0("<b>Method:</b> ", res$method, "<br/>",
           "<b>Samples:</b> ", paste(res$samples, collapse=", "), "<br/>",
           "<b>Conditions:</b> ", paste(res$conditions, collapse=", "), "<br/>",
           "<b>FDR cutoff:</b> ", res$fdr, "<br/>",
           "<b>LogFC cutoff:</b> ", res$lfc, "<br/>")
  })

  # normalized results (independent of test used)
  de.result.data <- function(res) {
    tab <- res$tab

    lfc <- switch(res$method,
                  "edgeR classic"=tab$logFC,
                  "edgeR_glm"=NULL,
                  "DESeq2"=tab$log2FoldChange)

    logCPM <- switch(res$method,
                  "edgeR classic"=tab$logCPM,
                  "edgeR_glm"=NULL,
                  "DESeq2"=log2(tab$baseMean))

    fdr <- switch(res$method,
                  "edgeR classic"=tab$FDR,
                  "edgeR_glm"=NULL,
                  "DESeq2"=tab$padj)

    data.frame(GeneID=rownames(tab), logFC=lfc, logCPM=logCPM, FDR=fdr, logFDR=-log10(fdr))
  }
  
  volcano.plot.fun <- function(res) {
    tab <- de.result.data(res)
    tab$Significant <- tab$FDR < res$fdr & abs(tab$logFC) >= res$lfc
    
    function() {
      ggplot(tab) +
        geom_point(aes(x=logFC, y=logFDR, col=Significant), size=0.75) +
        scale_color_manual(values=c("black", "red")) +
        geom_vline(xintercept=c(-res$lfc, res$lfc), lty="dashed") +
        theme(legend.position = "none")
    }
  }
  
  ma.plot.fun <- function(res) {
    tab <- de.result.data(res)
    tab$Significant <- tab$FDR < res$fdr & abs(tab$logFC) >= res$lfc

    function() {    
      ggplot(tab) + 
        geom_point(aes(x=logCPM, y=logFC, col=Significant), size=0.75) +
        scale_color_manual(values=c("black", "red")) +
        geom_vline(xintercept=0, lty="dashed") +
        theme(legend.position = "none")
    }
  }
  
  #------------------------------------------------------------------------------------------------------------  
  # aggs: add de.heat.plot.fun
  # it uses the rlog transformation of DESeq2 to plot DGE normalized and scaled 
  # or it uses the log2 CPM normalized and scaled counts of edgeR
  
  de.heat.plot.fun <- function(res) {
    tab <- de.result.data(res)
    tab$Significant <- tab$FDR < res$fdr & abs(tab$logFC) >= res$lfc
    mtx <- res$ge_mtx_trans
    de_genes <- tab[ tab$Significant %in% "TRUE", ] # subset the df for DGE
    topUp <- de_genes %>% 
      filter( logFC > 0 ) %>% 
      arrange( desc(logFC) ) %>% 
      pull( GeneID ) %>% 
      as.character(.)
    topDown <- de_genes %>% 
      filter( logFC < 0 ) %>% 
      arrange( logFC ) %>% 
      pull( GeneID ) %>% 
      as.character(.)
    if ( ( length(topUp) >= 25 ) & ( length(topDown) >= 25 ) ) {
      topUp <- topUp[1:25]
      topDown <- topDown[1:25]
    } else if ( ( length(topUp) >= 25 ) & ( length(topDown) < 25 ) ) {
      topUp <- topUp[1:25]
      topDown <- topDown
    } else if ( ( length(topUp) < 25 ) & ( length(topDown) >= 25 ) ) {
      topUp <- topUp
      topDown <- topDown[1:25]
    } else {
      topUp <- topUp
      topDown <- topDown
    }
    topDE <- c( topUp, topDown )
    #down_ge <- table(de_genes$logFC < 0 )["TRUE"]
    #top_dge <- table(de_genes$logFC > 0 )["TRUE"] 
    #colDataRow <- data.frame( "Regulation" = ifelse( de_genes$logFC < 0, paste0( "Downregulated (n=", down_ge, ")"), 
    #                                                 paste0("Upregulated (n=", top_dge, ")") ) ) # df for row annt
    colDataRow <- data.frame( "Regulation" = c( rep( "Upregulated", length(topUp) ), rep( "Downregulated", length(topDown) ) ) ) # df for row annt
    rownames(colDataRow) <- topDE
    mtx <- mtx[ topDE, ] # select only genes DGE
    mtx <- scale(mtx) # scale - z-scores 
    colDataCol <- data.frame( "Condition" = res$conditions ) # df for col annt
    rownames(colDataCol) <- colnames(mtx) # add row names to df
    heat_map_de <- heatmaply::ggheatmap(mtx , row_dend_left = TRUE,
                                        col_side_colors = colDataCol$Condition, 
                                        row_side_colors = colDataRow$Regulation, 
                                        showticklabels = c(TRUE, TRUE) )
    #heat_map_de <- pheatmap( mtx, 
    #                         show_rownames = FALSE, 
    #                         annotation_col = colDataCol, 
    #                         annotation_row = colDataRow )
    listOut <- list(
      mtx = mtx, 
      col_side_colors = colDataCol$Condition, 
      row_side_colors = colDataRow$Regulation,
      plot = heat_map_de
    )
    return(listOut)
  }
  #------------------------------------------------------------------------------------------------------------  
  
  #------------------------------------------------------------------------------------------------------------  
  # aggs: add dot.plot.top.up.fun and dot.plot.top.down.fun
  # it uses the normalized counts with the DESeq2 method
  # or it uses the CPM normalized of edgeR of the top down and up regulated genes
  
  dot.plot.top.up.fun <- function(res) {
    de_genes <- de.result.data(res)
    mtx <- res$ge_mtx_norm
    de_genes <- de_genes[ de_genes$FDR < res$fdr & abs(de_genes$logFC) >= res$lfc , ] # subset the df for DGE
    topUp <- de_genes %>% 
      filter( logFC > 0 ) %>% 
      arrange( desc(logFC) ) %>% 
      pull( GeneID ) %>% 
      as.character(.)
    topUp_len <- length( topUp )
    if ( topUp_len >= 10 ){
      dge_df <- mtx[ topUp[1:10], ] %>% 
        as.data.frame(.)
      conditions <- as.character( res$conditions )
      names( conditions ) <- colnames( dge_df )
      dge_df[ , "GeneID" ] <- factor( topUp[1:10], levels = topUp[1:10] ) 
      dge_df <- tidyr::gather( dge_df, "Samples", "Counts", colnames(dge_df)[1:(ncol(dge_df)-1)] )
      dge_df[ , "Condition" ] <- ifelse( dge_df$Samples %in% names( conditions[ conditions == "ConditionA" ] ), 
                                         "ConditionA", "ConditionB" )
      #browser()
      #readr::write_tsv(dge_df, "/home/agsousa/Desktop/data.frame.tsv")
      dot_plot <- 
        ggplot( dge_df, aes( x = Samples, y = Counts, color = Condition, fill = Condition ) ) + 
          geom_point() + 
          facet_wrap( ~ GeneID ) + 
          scale_y_log10() + ylab("Counts") + 
          #theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
          ggtitle("Top 10 upregulated genes")
    } else if ( (topUp_len > 0) & (topUp_len < 10) ) {
      dge_df <- mtx[ topUp, ] %>% 
        as.data.frame(.)
      conditions <- as.character( res$conditions )
      names( conditions ) <- colnames( dge_df )
      dge_df[ , "GeneID" ] <- topUp
      dge_df[ , "GeneID" ] <- factor( topUp, levels = topUp ) 
      dge_df <- tidyr::gather( dge_df, "Samples", "Counts", colnames(dge_df)[1:(ncol(dge_df)-1)] )
      dge_df[ , "Condition" ] <- ifelse( dge_df$Samples %in% names( conditions[ conditions == "ConditionA" ] ), 
                                         "ConditionA", "ConditionB" )
      dot_plot <- 
          ggplot( dge_df, aes( x = Samples, y = Counts, color = Condition, fill = Condition ) ) + 
          geom_point() + 
          facet_wrap( ~ GeneID ) + 
          scale_y_log10() + ylab("Counts") + 
          #theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
          ggtitle( paste0("Top ", topUp_len, " upregulated genes") )
    } else {
      dot_plot <- 
        ggplot() + 
        ggtitle("There were not found DGE upregulated to plot!")
    }
    
    return( dot_plot )
  }
  
  #--------------------------------------------------------------------------------------------------------------
  
  dot.plot.top.down.fun <- function(res) {
    de_genes <- de.result.data(res)
    mtx <- res$ge_mtx_norm
    de_genes <- de_genes[ de_genes$FDR < res$fdr & abs(de_genes$logFC) >= res$lfc , ] # subset the df for DGE
    topDown <- de_genes %>% 
      filter( logFC < 0 ) %>% 
      arrange( desc(logFC) ) %>% 
      pull( GeneID ) %>% 
      as.character(.)
    topDown_len <- length( topDown )
    if ( topDown_len >= 10 ){
      dge_df <- mtx[ topDown[1:10], ] %>% 
        as.data.frame(.)
      conditions <- as.character( res$conditions )
      names( conditions ) <- colnames( dge_df )
      dge_df[ , "GeneID" ] <- topDown[1:10]
      dge_df[ , "GeneID" ] <- factor( topDown[1:10], levels = topDown[1:10] ) 
      dge_df <- tidyr::gather( dge_df, "Samples", "Counts", colnames(dge_df)[1:(ncol(dge_df)-1)] )
      dge_df[ , "Condition" ] <- ifelse( dge_df$Samples %in% names( conditions[ conditions == "ConditionA" ] ), 
                                         "ConditionA", "ConditionB" )
      dot_plot <-
        ggplot( dge_df, aes( x = Samples, y = Counts, color = Condition, fill = Condition ) ) + 
          geom_point() + 
          facet_wrap( ~ GeneID ) + 
          scale_y_log10() + ylab("Counts") + 
          #theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
          ggtitle("Top 10 downregulated genes")
    } else if ( (topDown_len > 0) & (topDown_len < 10) ) {
      dge_df <- mtx[ topDown, ] %>% 
        as.data.frame(.)
      conditions <- as.character( res$conditions )
      names( conditions ) <- colnames( dge_df )
      dge_df[ , "GeneID" ] <- topDown
      dge_df[ , "GeneID" ] <- factor( topDown, levels = topDown ) 
      dge_df <- tidyr::gather( dge_df, "Samples", "Counts", colnames(dge_df)[1:(ncol(dge_df)-1)] )
      dge_df[ , "Condition" ] <- ifelse( dge_df$Samples %in% names( conditions[ conditions == "ConditionA" ] ), 
                                         "ConditionA", "ConditionB" )
      dot_plot <-
        ggplot( dge_df, aes( x = Samples, y = Counts, color = Condition, fill = Condition ) ) + 
          geom_point + 
          facet_wrap( ~ GeneID ) + 
          scale_y_log10() + ylab("Counts") + 
          #theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
          ggtitle( paste0("Top ", topDown_len, " downregulated genes") )
    } else {
      dot_plot <-
        ggplot() + 
        ggtitle("There were not found DGE downregulated to plot!")
    }
    return( dot_plot )
  }
  #------------------------------------------------------------------------------------------------------------ 
  
  
  #------------------------------------------------------------------------------------------------------------  
  # aggs: add pca_plot and heat_plot
  
  output$dim_plot <- renderPlot( {
    
    result()$plot_dim
    
  } )
  
  output$heat_plot <- renderPlot( {
    
    result()$plot_heat
    
  } )
  
  output$heat_de_plot <- renderPlot( {

    result()$de.heat.plot.fun$plot
    
  } )
  
  output$dot_dge_up_plot <- renderPlot( {

    result()$dot.plot.top.up.fun
    
  } )
  
  
  output$dot_dge_down_plot <- renderPlot( {
    
    result()$dot.plot.top.down.fun
    
  } )
  
  #------------------------------------------------------------------------------------------------------------
  
  # show the volcano plot
  output$volcano_plot <- renderPlot({
    result()$volcano.plot.fun()
  })
  
  # show the volcano plot
  output$plot_ma <- renderPlot({
    result()$ma.plot.fun()
  })
  
  # show the volcano plot
  output$plotly_plot <- renderPlotly({
    res <- result()
    tab <- de.result.data(res)
    
    tab$Significant <- tab$FDR < res$fdr & abs(tab$logFC) >= res$lfc
    
    if (input$radio_plotly_type == "Volcano plot") {
      plot_ly(tab, x=~logFC, y=~logFDR, color=~Significant, type="scatter", colors=c("black", "red"),
              text = ~paste("GeneID:", GeneID, "<br>logCPM:", logCPM, "<br>logFC:", logFC, '<br>FDR:', FDR))
    } else if (input$radio_plotly_type == "MA-plot") {
      plot_ly(tab, x=~logCPM, y=~logFC, color=~Significant, type="scatter", colors=c("black", "red"),
              text = ~paste("GeneID:", GeneID, "<br>logCPM:", logCPM, "<br>logFC:", logFC, '<br>FDR:', FDR))
    } else if (input$radio_plotly_type == "MDS/PCA") { # aggs
      plot_dim <- result()$plot_dim + 
        theme_bw()
      ggplotly( p = plot_dim )
    } else if (input$radio_plotly_type == "Heatmap - QC") { # aggs
      heatmaply::heatmaply( x = result()$corr_mtx, row_dend_left = TRUE, 
                 col_side_colors = result()$conditions)
    } else if (input$radio_plotly_type == "Heatmap - DE") { # aggs
      heat_de <-  result()$de.heat.plot.fun
      heatmaply::heatmaply( x = heat_de$mtx, row_dend_left = TRUE, 
                            col_side_colors = heat_de$col_side_colors, 
                            row_side_colors = heat_de$row_side_colors )
    } else if (input$radio_plotly_type == "Dotplot - up") { # aggs
      dot_plot_up <- result()$dot.plot.top.up.fun
      dot_plot_up$data[ dot_plot_up$data$Counts == 0, "Counts" ] <- 1e-06 # add this fake low no. 
      #to avoid infinite numbers when scaling to log10()
      dot_plot_up <- dot_plot_up + 
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      ggplotly( p = dot_plot_up )
    } else if (input$radio_plotly_type == "Dotplot - down") { # aggs
      dot_plot_down <- result()$dot.plot.top.down.fun
      dot_plot_down$data[ dot_plot_down$data$Counts == 0, "Counts" ] <- 1e-06 # add this fake low no. 
      #to avoid infinite numbers when scaling to log10()
      dot_plot_down <- dot_plot_down + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      ggplotly( p = dot_plot_down )
    } 
  })
  
  # do analysis when button is clicked
  observeEvent(input$btnGo, {
    res <- switch(input$sel_method,
                  "edger_classic"=edgeR_classic$result(),
                  "edgeR_glm"=NULL,
                  "deseq2"=deseq2$result())
    
    res$volcano.plot.fun <- volcano.plot.fun(res)
    res$ma.plot.fun <- ma.plot.fun(res)
    res$de.heat.plot.fun <- de.heat.plot.fun(res) # aggs
    res$dot.plot.top.up.fun <- dot.plot.top.up.fun(res) # aggs
    res$dot.plot.top.down.fun <- dot.plot.top.down.fun(res) # aggs
    
    if (is.na(res$error)) {
      results$data[[ res$id ]] <- res
    }
  })
  
  observeEvent(input$button_delete, {
    results$data[[ input$sel_result ]] <- NULL
  })

  # make the results panel show only if there are any results
  output$showResults <- reactive({
    # print("show results")
    # print(length(results$data) > 0)
    # 
    # length(results$data) > 0
    TRUE
  })
  outputOptions(output, "showResults", suspendWhenHidden = FALSE)
  
  # return results
  sessionData$de_results <- results
  
  return(sessionData)
}











