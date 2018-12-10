
source("mod_import_table.R")

#' UI function for table import module
tab_importUI <- function(id, datasets) {
  ns <- NS(id)
  
  panels.ui <- tagList(h4("Summary"),
                       htmlOutput(ns("table_summary1")),
                       hr(),
                       h4("Preview"),
                       tableOutput(ns("table_preview")),
                       hr(),
                       h4("Sample Metadata"),
                       tableOutput(ns("table_metadata")))
  
  import.ui <- tagList(
    radioButtons(ns("radioImport"), label = "Import from",
                 choices = list("File upload" = "file", "Built-in dataset" = "dataset"), 
                 selected = "dataset"),
    tags$hr(),
    conditionalPanel(condition = paste0("input['", ns("radioImport"), "']", " == 'file'"), 
                     tagList(
                       h4("Count Data"),
                       helpText("Gene names should be on the first column,",
                                "and every other column should contain the raw count data for each sample."),
                       mod_import_tableUI(ns("table_import")),
                       hr(),
                       h4("Sample Metadata (Optional)"),
                       helpText("Sample names should appear in the first column."),
                       mod_import_tableUI(ns("table_metadata_import"))
                     )),
    conditionalPanel(condition = paste0("input['", ns("radioImport"), "']", " == 'dataset'"), 
                     mod_import_datasetUI(ns("dataset_import"), datasets))
  )
  
  sidebarLayout(
    sidebarPanel(import.ui),
    mainPanel(panels.ui))
}

#' Server function for table loader module
#' 
#' @return A dataframe as a reactive value.
tab_importServer <- function(input, output, session, datasets, sessionData) {
  
  fileImportData <- callModule(mod_import_tableServer, "table_import", stringsAsFactors = FALSE)
  metadataImportData <- callModule(mod_import_tableServer, "table_metadata_import", stringsAsFactors = FALSE)
  
  datasetImportData <- callModule(mod_import_datasetServer, "dataset_import", datasets)
  
  dataframe <- reactive({
    req(input$radioImport)
    
    if (input$radioImport == "file") {
      req(fileImportData$dataframe())
    } else {
      req(datasetImportData$dataframe())
    }
    
    df <- switch(input$radioImport,
           "file"=fileImportData$dataframe(),
           "dataset"=datasetImportData$dataframe())
    
    rownames(df) <- df[, 1]
    df[, 1] <- NULL
    
    return(df)
  })

  metadata <- reactive({
    req(input$radioImport)
    
    df <- switch(input$radioImport,
           "file"=metadataImportData$dataframe(),
           "dataset"=datasetImportData$metadata())
    
    rownames(df) <- df[, 1]
    df[, 1] <- NULL
    
    return(df)
  })
  
  #observeEvent(dataframe())
  
  # show the summary
  output$table_summary1 <- renderUI({
    df <- dataframe()
    
    txt <- paste(paste("Number of columns:", ncol(df)),
                 paste("Number of rows:", nrow(df)),
                 paste("Number of numeric columns:", sum(apply(df, 1, is.numeric))),
                 paste("Number of text columns:", sum(!apply(df, 1, is.numeric))),
                 paste("Columns with NA's:", sum(apply(df, 2, function(x) sum(is.na(x)) > 0))),
                 paste("Rows with NA's:", sum(apply(df, 1, function(x) sum(is.na(x)) > 0))),
                 sep="<br/>")
    
    return (HTML(txt))
  })

  # show the table
  output$table_preview <- renderTable({
    df <- dataframe()
    
    validate(need(df, "Please import a table."))
    
    nc <- min(ncol(df), 6)
    nr <- min(nrow(df), 6)
    
    return(df[ 1:nr, 1:nc ])
  }, rownames = TRUE)
  
  
  output$table_metadata <- renderTable({
    df <- metadata()
    
    validate(need(df, "Please import a table."))
    
    return(df)
  }, rownames = TRUE)
  

  
  sessionData$dataframe <- dataframe
  sessionData$metadata <- metadata
  
  return(sessionData)
}


