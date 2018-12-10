







#' UI function for the table filter module.
tableRenameUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    uiOutput(ns("rename_ui"))
  )
}

#' Server function for the table filter module.
tableRenameServer <- function(input, output, session, dataframe) {
  
  
  output$rename_ui <- renderUI({
    ns <- session$ns
    cnames <- colnames(dataframe())
    
    tagList(
      lapply(cnames, function(x) {
        textInput(ns(paste0("input-", x)), label = x, value = x)
      })
    )
  })
  
  oldnames <- reactive({
    colnames(dataframe())
  })
  
  newnames <- reactive({
    res <- sapply(oldnames(), function(x) {
      input[[ paste0("input-", x) ]]
    })
    
    w <- which(res == "NULL")
    res[ w ] <- oldnames()[w]
    
    print(w)
    print(res)
    
    return(res)
  })
  
  renamed <- reactive({
    df <- dataframe()
    
    colnames(df) <- newnames()
    
    return (df)
  })
  
  # return the reactive value
  return (renamed)
}














#' UI function for the table filter module.
mod_filter_tableUI <- function(id, dataframe=NULL) {
  ns <- NS(id)
  
  if (is.null(dataframe)) {
    choices <- NULL
  } else {
    choices <- colnames(dataframe)
  }
  
  col.filter <- tagList(
    selectizeInput(ns("columns"), "Use columns", choices = choices, selected = choices, multiple = TRUE, width="100%"),
    actionButton(ns("btnAll"), "Select All"),
    actionButton(ns("btnNone"), "Select None")
    
  )
  
  row.filter <- tagList(
    checkboxInput(ns("chkSample"), label = "Sample rows", value = FALSE),
    conditionalPanel("input.chkSample == true", ns = ns,
                     numericInput(ns("numSampleSize"), label = "Sample size", value = 1000))
  )
  
  # tagList(
  #   useShinyjs(),
  #   bsCollapse(id=ns("type"), multiple=TRUE, open = c("Columns", "Rows"),
  #              bsCollapsePanel("Columns", col.filter, style = "default"),
  #              bsCollapsePanel("Rows", row.filter, style = "default")
  #   )
  # )
  
  tagList(
    col.filter,
    row.filter
    )

}

#' Server function for the table filter module.
mod_filter_tableServer <- function(input, output, session, dataframe) {
  
  # select all columns if raw dataframe changes
  observe({
    dataframe()
    #input$columns
    
    req(dataframe())
    
    print("update selectize input choices")
    
    valid.choices <- colnames(dataframe())
    
    print(valid.choices)
    
    updateSelectizeInput(session, "columns", choices=valid.choices, selected=valid.choices, server=TRUE)
  })
  
  # # update column chooser if numeric filter changes
  # observe({
  #   print("update selectize input on numeric change")
  #   input$numeric # observe this
  #   
  #   # store selected columns (isolating the column selector)
  #   selected <- isolate(input$columns)
  #   valid.choices <- isolate(choices())
  #   
  #   if (length(valid.choices) > 0)
  #     selected <- intersect(selected, valid.choices)
  #   
  #   updateSelectizeInput(session, "columns", choices=valid.choices, selected=selected, server=TRUE)
  # })
  
  # select all columns
  observeEvent(input$btnAll, {
    print("select all")
    valid.choices <- colnames(dataframe())
    
    updateSelectizeInput(session, "columns", choices=valid.choices, selected=valid.choices, server=TRUE)
  })
  
  # remove all columns
  observeEvent(input$btnNone, {
    print("select none")
    valid.choices <- colnames(dataframe())
    
    updateSelectizeInput(session, "columns", choices=valid.choices, selected=NULL, server=TRUE)
  })
  
  # get the filtered data.frame
  filtered <- reactive({
    # dependencies must come before isolate... why???
    input$columns
    input$chkSample
    input$numSampleSize
    
    df <- isolate(dataframe())
    
    req(df)
    #validate(need(length(input$columns) > 0, "Error: At least one column must be selected."))
    validate(need(!is.na(input$numSampleSize), "Error: Sample size must be numeric."))
    
    print("filter dataframe")
    
    # select columns
    filtered.df <- df[ , match(input$columns, colnames(df)), drop = FALSE ]
    
    # sample rows
    if (input$chkSample == TRUE) {
      set.seed(42)
      
      sample.size <- min(nrow(filtered.df), input$numSampleSize)
      w <- sample(1:nrow(filtered.df), sample.size)
      
      filtered.df <- filtered.df[ w,, drop = FALSE ]
    }
    
    return (filtered.df)
  })
  
  # return the reactive value
  return (filtered)
}

