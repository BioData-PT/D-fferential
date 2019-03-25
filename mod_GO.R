library(png)

mod_GOUI <- function(id, ensembl) {
  ns <- NS(id)
  
  ensembl.choices <- sort(listDatasets(ensembl)$description)
  
  go.options <- tagList(
    h4("GO Enrichment Options"),
    selectInput(ns("sel_annotation"), label = "Annotation", choices = ensembl.choices),
    selectInput(ns("sel_attribute"), label = "Gene IDs", choices = NULL),
    selectInput(ns("sel_de_result"), label = "Geneset Source", choices = NULL),
    radioButtons(ns("radio_genes"), "Population", 
                 choices = list("Up-regulated"="up", "Down-regulated"="down", "Differentially expressed"="all")),
    hr(),
    textInput(ns("text_name"), "Name for this analysis"),
    actionButton(ns("btn_go"), "Go!", width="100%"))
  
  go.panel <- tabsetPanel(
    tabPanel("Annotation", DTOutput(ns("tab_annotation"))),
    #tabPanel("Universe", verbatimTextOutput(ns("text_universe")))
    tabPanel("Biolocal Process", 
             uiOutput(ns("ui_bp_graph")),
             hr(),
             DTOutput(ns("tab_bp"))),
    tabPanel("Cellular Components", 
             uiOutput(ns("ui_cc_graph")),
             hr(),
             DTOutput(ns("tab_cc"))),
    tabPanel("Molecular Function", 
             uiOutput(ns("ui_mf_graph")),
             hr(),
             DTOutput(ns("tab_mf")))
  )
  
  main.panel.ui <- tagList(
    selectInput(ns("selResult"), label = "Show Result", choices = NULL),
    htmlOutput(ns("html_res_params")),
    hr(),
    go.panel,
    verbatimTextOutput(ns("debug"))
  )
  
  # sidebarLayout(
  #   sidebarPanel(go.options),
  #   mainPanel(main.panel.ui)
  # )
  fluidRow(
    box(go.options, width=4),
    box(main.panel.ui, width=8)
  )
  
}

mod_GOServer <- function(input, output, session, sessionData, ensembl) {
  
  print("starting mod go")
  ensembl_datasets <- listDatasets(ensembl)
  print(head(ensembl_datasets))
  
  results <- reactiveValues()
  
  result <- reactive({
    req(input$selResult)
    
    results[[ input$selResult ]]
  })
  
  # start with the Go button disabled
  shinyjs::disable("btn_go")

  annotation_filename <- reactive({
    req(input$sel_annotation, input$sel_attribute)
    
    annotation_description <- isolate(input$sel_annotation)
    attribute <- isolate(input$sel_attribute)
    
    dataset.name <- ensembl_datasets$dataset[ match(annotation_description, ensembl_datasets$description) ]
    dataset.pathname <- file.path("data_persist", paste0(dataset.name, "_", attribute, ".txt"))
    
    if (!file.exists(dataset.pathname)) {
      withProgress(message = paste0('Downloading GO annotation for ', annotation_description, "...") , {
        mart <- useEnsembl(biomart="ensembl", dataset=dataset.name)
        bm <- getBM(attributes = c(attribute, "go_id"), mart = mart)
        bm <- bm[ which(bm$go_id != ""), ]
        
        write.table(bm, dataset.pathname, quote = FALSE, row.names = FALSE, sep="\t")
      })
    } 

    return(dataset.pathname)
  })
  
  ensembl_dataset <- reactive({
    dataset <- read.table(annotation_filename())
    
    return(dataset)
  })
  
  observe({
    updateSelectInput(session, "sel_de_result", choices=names(sessionData$de_results$data))
  })
  
  # Observe a change in the selected annotation and update available attributes
  observe({
    req(input$sel_annotation)
    
    dataset.name <- ensembl_datasets$dataset[ match(input$sel_annotation, ensembl_datasets$description) ]
    mart <- useEnsembl(biomart="ensembl", dataset=dataset.name)
        
    atts <- listAttributes(mart = mart)
    atts.gene <- grep("gene", atts$name, value=TRUE)
    
    possible.choices <- c("ensembl_gene_id", "entrezgene")
    
    atts.gene <- intersect(atts.gene, possible.choices)
    
    updateSelectInput(session, "sel_attribute", choices=atts.gene)
  })
  
  # enable the Go button when ready
  observe({
    if (input$sel_de_result == "")
      shinyjs::disable("btn_go")
    else {
      shinyjs::enable("btn_go")
    }
  })

  output$tab_annotation <- renderDT({
    #ensembl_dataset()
    result()$annotation
  })
  
  output$text_universe <- renderText({
    paste(result()$universe, collapse="\n")
  })
  
  obo.file <- reactive({
    obo.filename <- "data_persist/go.obo"
    
    if (!file.exists(obo.filename)) {
      withProgress(message = 'Downloading GO Core ontology...', {
        download.file("http://purl.obolibrary.org/obo/go.obo", obo.filename)
      })
    }
    
    if (file.exists(obo.filename))
      return(obo.filename)
    else
      return(NULL)
  })
  
  de_result <- reactive({
    req(input$sel_de_result)
    
    sessionData$de_results$data[[ input$sel_de_result ]]
  })
  
  # update radio buttons when de_result is changed
  observe({
    de <- de_result()

    if (!is.null(de)) {
      choices <- list()
      choices[[ paste0("Up-regulated (", length(de$genes.up), ")") ]] <- "up"
      choices[[ paste0("Down-regulated (", length(de$genes.down), ")") ]] <- "down"
      choices[[ paste0("Differentially expressed (", length(de$genes.diff), ")") ]] <- "all"
    } else {
      choices = list("Up-regulated"="up", "Down-regulated"="down", "Differentially expressed"="all")
    }
    
    print(choices)
    
    updateRadioButtons(session, "radio_genes", choices = choices)
  })
  
  observeEvent(input$btn_go, {
    print ("Running GOEnrichment")
    req(obo.file(), annotation_filename())
    
    # universe
    universe <- rownames(de_result()$table)
    universe.filename <- tempfile()
    writeLines(universe, universe.filename)
    
    # population
    study <- switch(input$radio_genes,
                    "up" = de_result()$genes.up,
                    "down" = de_result()$genes.down,
                    "all" = de_result()$genes.diff)
    
    study.filename <- tempfile()
    writeLines(study, study.filename)
    
    mf.filename <- tempfile(pattern="GO_MF")
    cc.filename <- tempfile(pattern="GO_CC")
    bp.filename <- tempfile(pattern="GO_BP")
    
    mf.graph.filename <- tempfile(pattern="GO_GRAPH_MF")
    cc.graph.filename <- tempfile(pattern="GO_GRAPH_CC")
    bp.graph.filename <- tempfile(pattern="GO_GRAPH_BP")
    
    CMD <- paste("java -jar bin/GOEnrichment.jar",
                 "-g", obo.file(),
                 "-a", annotation_filename(),
                 "-e",
                 "-so",
                 "-c Benjamini-Hochberg",
                 "-p", universe.filename,
                 "-s", study.filename, 
                 "-mfr", mf.filename,
                 "-ccr", cc.filename,
                 "-bpr", bp.filename,
                 "-mfg", mf.graph.filename,
                 "-ccg", cc.graph.filename,
                 "-bpg", bp.graph.filename)
    
    withProgress(message = "Running GOEnrichment...", {
      system(CMD)
    })
    
    n <- length(names(results)) + 1
    
    if (nchar(input$text_name) > 0) {
      id <- input$text_name
    } else {
      id <- paste("GO enrichment", n)
    }
    
    result <- list(
      id = id,
      annotation = ensembl_dataset(),
      annotation.name = input$sel_annotation,
      attribute.name = input$sel_attribute,
      geneset.source = input$sel_de_result,
      geneset.type = input$radio_genes,
      universe = universe,
      study = study,
      mf.graph.filename = mf.graph.filename,
      cc.graph.filename = cc.graph.filename,
      bp.graph.filename = bp.graph.filename,
      mf.graph = openssl::base64_encode(readBin(mf.graph.filename, "raw", file.info(mf.graph.filename)[1, "size"])),
      cc.graph = openssl::base64_encode(readBin(cc.graph.filename, "raw", file.info(cc.graph.filename)[1, "size"])),
      bp.graph = openssl::base64_encode(readBin(bp.graph.filename, "raw", file.info(bp.graph.filename)[1, "size"])),
      mf.tab = read.table(mf.filename, sep="\t", header=TRUE, comment.char = "", quote = ""),
      cc.tab = read.table(cc.filename, sep="\t", header=TRUE, comment.char = "", quote = ""),
      bp.tab = read.table(bp.filename, sep="\t", header=TRUE, comment.char = "", quote = "")
    )
    
    results[[ id ]] <- result
  })
  
  # observe new results and update the select input and name for next analysis  
  observe({
    updateTextInput(session, "text_name", value = paste("GO enrichment", length(names(results)) + 1))
    updateSelectInput(session, "selResult", choices = names(results), selected = rev(names(results))[1])
  })
  
  output$html_res_params <- renderText({
    res <- result()
    
    paste0("<b>Annotation:</b> ", res$annotation.name, "<br/>",
           "<b>GeneID Format:</b> ", res$attribute.name, "<br/>",
           "<b>Geneset source:</b> ", res$geneset.source, "<br/>",
           "<b>Background Genes:</b> ", length(res$universe), "<br/>",
           "<b>Population Genes:</b> ", length(res$study), "<br/>")
  })

  image_info <- function(filename) {
    list(src = filename, contentType = 'image/png')
  }
  
  output$ui_bp_graph <- renderUI({
    img <- readPNG(result()$bp.graph.filename)
    imageOutput(session$ns("img_bp_graph"), width=dim(img)[2], height=dim(img)[1])
  })
  output$ui_cc_graph <- renderUI({
    img <- readPNG(result()$cc.graph.filename)
    imageOutput(session$ns("img_cc_graph"), width=dim(img)[2], height=dim(img)[1])
  })
  output$ui_mf_graph <- renderUI({
    img <- readPNG(result()$mf.graph.filename)
    imageOutput(session$ns("img_mf_graph"), width=dim(img)[2], height=dim(img)[1])
  })
  
  output$img_bp_graph <- renderImage(image_info(result()$bp.graph.filename), deleteFile = FALSE)
  output$img_cc_graph <- renderImage(image_info(result()$cc.graph.filename), deleteFile = FALSE)
  output$img_mf_graph <- renderImage(image_info(result()$mf.graph.filename), deleteFile = FALSE)
  
  output$tab_cc <- renderDT(result()$cc.tab)
  output$tab_bp <- renderDT(result()$bp.tab)
  output$tab_mf <- renderDT(result()$mf.tab)
  
  sessionData$go_results <- results
  
  return(sessionData)
}

