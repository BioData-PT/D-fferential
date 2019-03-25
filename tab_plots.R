
source("mod_filter_table.R")

# # Density plot UI
# #################################################
# mod_plot_densityUI <- function(id, options) {
#   ns <- NS(id)
#   
#   col.choices <- c("None", "SampleID")
#   if (!is.null(options$metadata)) {
#     col.choices <- c(col.choices, colnames(options$metadata))
#   }
#   
#   tagList(
#     bsCollapsePanel("Plot options", 
#       inputPanel(mod_filter_tableUI(ns("mod_filter"), options$dataframe)),
#       inputPanel(
#         selectInput(ns("sel_scale"), label = "Scale", 
#                     choices = c("Linear"="identity", 
#                                 "Log2"="log2", 
#                                 "Log10"="log10"), selected = "identity"),
#         selectInput(ns("sel_col"), "Color", choices = col.choices),
#         selectInput(ns("sel_lty"), "Line type", choices = col.choices),
#         checkboxInput(ns("chk_legend"), label = "Show legend", value = TRUE)
#       )
#     ),
#     plotOutput(ns("my_plot"))
#   )
# }
# 
# # Density plot server
# #################################################
# mod_plot_densityServer <- function(input, output, session, options) {
#   # this makes a copy of the provided data
#   dataframe <- options$dataframe
#   metadata <- options$metadata
#   
#   filtered <- callModule(mod_filter_tableServer, "mod_filter", reactive(dataframe))
# 
#   data.points <- reactive({
#     req(filtered())
#     
#     validate(need(ncol(filtered()) > 0, "No columns."))
#     
#     mdf <- melt(filtered())
#     colnames(mdf) <- c("SampleID", "value")
#     
#     if (!is.null(metadata)) {
#       mdf <- cbind(mdf, metadata[ match(mdf$SampleID, rownames(metadata)), ])
#     }
#     
#     return(mdf)
#   })
#   
#   density.plot <- reactive({
#     mdf <- data.points()
#     
#     color_by <- if (input$sel_col == "None") NULL else input$sel_col
#     lty_by <- if (input$sel_lty == "None") NULL else input$sel_lty
#     trans <- input$sel_scale
#     legend <- input$chk_legend
# 
#     plot.fun <- function() {
#       p <- ggplot(mdf, aes_string(x="value", col=color_by, lty=lty_by)) + 
#         geom_density() +
#         scale_x_continuous(trans=trans)
#       
#       if (legend == FALSE) {
#         p <- p + theme(legend.position = "none")
#       }
#       
#       return (p)
#     }
#     
#     return(plot.fun)
#   })
# 
#   output$my_plot <- renderPlot({
#     density.plot()()
#   })
#   
#   return(list(
#     plot.fun = density.plot
#   ))
# }

# Boxplot plot UI
#################################################
mod_plot_boxplotUI <- function(id, options) {
  ns <- NS(id)
  
  col.choices <- c("None", "SampleID")
  if (!is.null(options$metadata)) {
    col.choices <- c(col.choices, colnames(options$metadata))
  }
  
  tagList(
    bsCollapsePanel("Options", 
      inputPanel(mod_filter_tableUI(ns("mod_filter"), options$dataframe)),
      inputPanel(
        selectInput(ns("sel_scale"), label = "Scale", 
                    choices = c("Linear"="identity", 
                                "Log2"="log2", 
                                "Log10"="log10"), selected = "identity"),
        selectInput(ns("sel_group"), "Group by", choices = col.choices[ -1 ]),
        selectInput(ns("sel_col"), "Color", choices = col.choices),
        selectInput(ns("sel_lty"), "Line type", choices = col.choices),
        checkboxInput(ns("check_notch"), "Show notch"),
        checkboxInput(ns("chk_legend"), label = "Show legend", value = TRUE)
      )
    ),
    plotOutput(ns("my_plot"))
  )
}

# Boxplot plot server
#################################################
mod_plot_boxplotServer <- function(input, output, session, options) {
  # this makes a copy of the provided data
  dataframe <- options$dataframe
  metadata <- options$metadata
  
  filtered <- callModule(mod_filter_tableServer, "mod_filter", reactive(dataframe))
  
  data.points <- reactive({
    req(filtered())
    
    validate(need(ncol(filtered()) > 0, "No columns."))
    
    mdf <- melt(filtered())
    colnames(mdf) <- c("SampleID", "value")
    
    if (!is.null(metadata)) {
      mdf <- cbind(mdf, metadata[ match(mdf$SampleID, rownames(metadata)), ])
    }
    
    return(mdf)
  })
  
  boxplot.plot <- reactive({
    mdf <- data.points()
    
    group_by <- input$sel_group
    color_by <- if (input$sel_col == "None") NULL else input$sel_col
    lty_by <- if (input$sel_lty == "None") NULL else input$sel_lty
    trans <- input$sel_scale
    legend <- input$chk_legend
    
    plot.fun <- function() {
      p <- ggplot(mdf, aes_string(y="value", x=group_by, fill=color_by, lty=lty_by)) + 
        geom_boxplot() +
        scale_y_continuous(trans=trans)
      
      if (legend == FALSE) {
        p <- p + theme(legend.position = "none")
      }
      
      return (p)
    }
    
    return(plot.fun)
  })
  
  output$my_plot <- renderPlot({
    boxplot.plot()()
  })
  
  return(list(
    plot.fun = boxplot.plot
  ))
}








# mod_column_chooserUI <- function(id, datasets) {
#   ns <- NS(id)
# 
#   tagList(
#     selectInput(ns("sel_dataset"), "Dataset", choices = names(datasets)),
#     selectInput(ns("sel_column"), "Variable", choices = colnames(datasets[[1]]))
#   )
# }
# 
# mod_column_chooserServer <- function(input, output, session, datasets) {
#   data <- reactive({
#     dataset <- datasets[[ input$sel_dataset ]]
#     
#     return(dataset[ , input$sel_column ])
#   })
#   
#   return(data)
# }
# 
# # Scatterplot plot UI
# #################################################
# mod_plot_scatterplotUI <- function(id, options) {
#   ns <- NS(id)
#   
#   # col.choices <- c("None", "SampleID")
#   # if (!is.null(options$metadata)) {
#   #   col.choices <- c(col.choices, colnames(options$metadata))
#   # }
#   
#   datasets <- options$datasets
#   
#   
#   tagList(
#     bsCollapsePanel("Options", 
#                     inputPanel(mod_column_chooserUI(ns("mod_columns"), datasets))
#                     #inputPanel(mod_filter_tableUI(ns("mod_filter"), options$dataframe))
#                     # inputPanel(
#                     #   selectInput(ns("sel_scale"), label = "Scale", 
#                     #               choices = c("Linear"="identity", 
#                     #                           "Log2"="log2", 
#                     #                           "Log10"="log10"), selected = "identity"),
#                     #   selectInput(ns("sel_group"), "Group by", choices = col.choices[ -1 ]),
#                     #   selectInput(ns("sel_col"), "Color", choices = col.choices),
#                     #   selectInput(ns("sel_lty"), "Line type", choices = col.choices),
#                     #   checkboxInput(ns("check_notch"), "Show notch"),
#                     #   checkboxInput(ns("chk_legend"), label = "Show legend", value = TRUE)
#                     # )
#     ),
#     plotOutput(ns("my_plot"))
#   )
# }

# Scatterplot plot server
#################################################
mod_plot_scatterplotServer <- function(input, output, session, options) {

  
  output$my_plot <- renderPlot({
    #boxplot.plot()()
    NULL
  })
  
  # return(list(
  #   plot.fun = boxplot.plot
  # ))
}






# Plots tab UI
#################################################
tab_plotsUI <- function(id) {
  ns <- NS(id)
 
  main.panel.ui <- tagList(
    div(id = "elementsHead")
  )
  
  sidebarLayout(
    sidebarPanel(
      selectInput(ns("sel_plot_type"), label = "Type", 
                  choices = c("Scatterplot", "Boxplot", "Density")),
      actionButton(ns("btn_add"), label = "Add plot")),
    mainPanel(main.panel.ui))
}

# Plots tab Server
#################################################
tab_plotsServer <- function(input, output, session, sessionData) {
  plot.creation.data <- reactive({
    datasets <- list("Counts" = sessionData$dataframe())
    
    list(
      "Scatterplot" = list(
        ui.fun = mod_plot_scatterplotUI,
        server.fun = mod_plot_scatterplotServer,
        options = list(
          datasets = datasets,
          metadata = sessionData$metadata(),
          type = "Scatterplot"
        )),
      "Density" = list(
        ui.fun = mod_plot_densityUI,
        server.fun = mod_plot_densityServer,
        options = list(
          dataframe = sessionData$dataframe(),
          metadata = sessionData$metadata(),
          type = "Density"
        )),
      "Boxplot" = list(
        ui.fun = mod_plot_boxplotUI,
        server.fun = mod_plot_boxplotServer,
        options = list(
          dataframe = sessionData$dataframe(),
          metadata = sessionData$metadata(),
          type = "Boxplot"
        ))
    )
  })
  
  elem_list <- reactiveValues()
  remove_observers <- list()
  
  insert.plot <- function(id, title, plot.options) {
    remove.btn.id <- paste0("btn_remove_", id)
    
    plot.ui <- plot.options$ui.fun(session$ns(id), plot.options$options)
    plot.ui <- bsCollapse(id = session$ns(paste0(id, "_collapse")), open = title, 
                          bsCollapsePanel(title = title, 
                                          actionButton(session$ns(remove.btn.id), label="Remove plot"),
                                          plot.ui))

    insertUI(selector = "#elementsHead",
             ui = plot.ui)
    
    remove_observers[[ remove.btn.id ]] <- observeEvent(input[[ remove.btn.id ]], {
      remove_observers[[ remove.btn.id ]] <- NULL
      
      removeUI(selector = paste0('#', session$ns(paste0(id, "_collapse"))))
    })

    callModule(plot.options$server.fun, id = id, plot.options$options)
  }
  
  observeEvent(input$btn_add, {
    n <- length(names(elem_list)) + 1
    id <- paste0("elem_", n)

    plot.options <- plot.creation.data()[[ input$sel_plot_type ]]
    
    elem_list[[ id ]] <- insert.plot(id, input$sel_plot_type, plot.options)
  })
  
  sessionData$plots <- elem_list
  
  return(sessionData)
}




