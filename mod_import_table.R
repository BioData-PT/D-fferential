

#' UI function for table import module
mod_import_datasetUI <- function(id, datasets) {
  ns <- NS(id)
  
  tagList(
    selectInput(ns("selDataset"), label = "Dataset", choices = names(datasets)),
    tags$b("Description"),
    textOutput(ns("txtDescription"))
  )
}

#' Server function for table loader module
#' 
#' @return A dataframe as a reactive value.
mod_import_datasetServer <- function(input, output, session, datasets) {
  output$txtDescription <- renderText({ 
    print(as.character(datasets[[ input$selDataset ]]$description))
  })
  
  dataframe <- reactive({
    datasets[[ input$selDataset ]]$dataframe
  })
  
  metadata <- reactive({
    datasets[[ input$selDataset ]]$metadata
  })
  
  name <- reactive({
    input$selDataset
  })
  
  list(dataframe=dataframe,
       metadata=metadata,
       name=name)
}




#' UI function for table import module
mod_import_tableUI <- function(id) {
  ns <- NS(id)
  
  tagList(
  #  wellPanel(
    fileInput(ns("file"), "", width = "100%"),
    checkboxInput(ns("heading"), "Has heading", value = TRUE),
    fluidRow(
        column(6, selectInput(ns("sep"), 
                              "Separator", 
                              c("Space" = " ",
                                "Tab" = "\t",
                                "Comma" = ",",
                                "Semicolon" = ";"), 
                              selected = "\t")),
        column(6, selectInput(ns("quote"), 
                              "Quote", 
                              c("None" = "",
                                "Double quote" = "\"",
                                "Single quote" = "'"), 
                              selected="None"))
      )
    )
}

#' Server function for table loader module
#' 
#' @return A dataframe as a reactive value.
mod_import_tableServer <- function(input, output, session, stringsAsFactors) {
  # get the file
  fileHandle <- reactive({
    validate(need(input$file, message = FALSE))
    
    input$file
  })
  
  # parse into a data.frame
  dataframe <- reactive({
    read.table(fileHandle()$datapath,
             header = input$heading, 
             sep = input$sep,
             quote = input$quote,
             stringsAsFactors = stringsAsFactors)
  })
  
  name <- reactive({
    fileHandle()$filename
  })
  
  return(list(dataframe=dataframe,
              name=name))
}
