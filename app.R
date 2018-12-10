# Differential expression anlysis tool

# devtools::install_github('rstudio/DT')

library(DT)

library(shiny)
library(shinyBS)
library(shinyjs)

library(gridExtra)

library(biomaRt)
library(edgeR)
library(DESeq2)


# options(shiny.reactlog=TRUE)
# options(shiny.trace=TRUE)

#source("tableFilterModule.R")
#source("tableTransformModule.R")

source("tab_import.R")
source("tab_plots.R")
source("mod_filter_table.R")
source("tab_de.R")
source("mod_GO.R")
source("tab_export.R")

options(shiny.maxRequestSize=1*1024^3)

# example datasets
example.datasets <- list(
  trapnell = list(
    dataframe = read.table("data/trapnell_counts.tab", header=TRUE),
    description = "Trapnell in-silico dataset."),
  lactate = list(
    dataframe = read.table("data/edgeR_example4_GSE60450_Lactation-GenewiseCounts.tab", header=TRUE)[, -2],
    description = "Lactogenesis.",
    metadata = read.table("data/edgeR_example4_GSE60450_Lactation_metadata.tab", header=TRUE)
    ),
  tuch = list(
    dataframe = (function() { 
      ids <- c("N8", "N33", "N51", "T8", "T33", "T51")
      filenames <- file.path("data", paste0("Tuch_", ids, ".tab"))
      tab <- do.call(cbind, lapply(filenames, function(x) read.table(x)[, 2]))
      
      colnames(tab) <- ids
      
      tab <- cbind(ID=read.table(filenames[1])[, 1], tab)
      
      return(as.data.frame(tab))
    })(),
    description = "Tuch et al."
  )
)

ensembl <- useEnsembl(biomart="ensembl")

ui <- navbarPage(id = "tabs",
                 theme = "bootswatch-3.3.1/simplex/bootstrap.css",
                 title = "D-fferrential", 
                 tabPanel("Import", tab_importUI("tab_import", example.datasets)),
                 tabPanel("Plots", tab_plotsUI("tab_plots")),
                 tabPanel("Differential Expression", tab_deUI("tab_de")),
                 tabPanel("GO Enrichment", mod_GOUI("mod_GO", ensembl)),
                 tabPanel("Export Analysis", tab_exportUI("tab_export"))
)

#' Main application server function 
server <- function(input, output, session) {
  # hideTab("tabs", target = "Filter")
  
  sessionData <- list(
    dataframe = NULL
  )
  
  # import tab
  sessionData <- callModule(tab_importServer, "tab_import", example.datasets, sessionData)

  # plots tab
  sessionData <- callModule(tab_plotsServer, "tab_plots", sessionData)
  
  # DE tab
  sessionData <- callModule(tab_deServer, "tab_de", sessionData)
  
  # GO tab
  sessionData <- callModule(mod_GOServer, "mod_GO", sessionData, ensembl)
  
  # export
  sessionData <- callModule(tab_exportServer, "tab_export", sessionData)

  # observeEvent(sessionData$dataframe(), {
  #   req(sessionData$dataframe())
  #   
  #   showTab("tabs", target = "Filter")
  # })
  
}

shinyApp(ui, server)

