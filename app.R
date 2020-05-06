# Differential expression anlysis tool

# devtools::install_github('rstudio/DT')

library(DT)

library(shiny)
library(shinyBS)
library(shinyjs)
library(shinydashboard)

library(gridExtra)

library(biomaRt)
library(edgeR)
library(DESeq2)

#library(ggvis)
library(plotly)
#library("heatmaply") # aggs: it requires this package
#library("ggrepel) # aggs: it requires this package for label annotation
library(pheatmap)

# options(shiny.reactlog=TRUE)
# options(shiny.trace=TRUE)

#source("tableFilterModule.R")
#source("tableTransformModule.R")

library(ggplot2)
library(reshape2)

#source("mod_filter_table.R")
source("tab_import.R")
#source("tab_visualize.R") # aggs: create a new one 
#source("tab_plots.R") # aggs: change 'tab_plots.R' --> 'tab_visualize.R  '
source("tab_de.R")
source("mod_GO.R")
source("tab_export.R")

options(shiny.maxRequestSize=1*1024^3)

# # example datasets
# example.datasets <- list(
#   trapnell = list(
#     dataframe = read.table("data/trapnell_counts.tab", header=TRUE, row.names = 1),
#     description = "Trapnell in-silico dataset."),
#   lactate = list(
#     dataframe = read.table("data/edgeR_example4_GSE60450_Lactation-GenewiseCounts.tab", header=TRUE, row.names = 1)[, -1],
#     description = "Lactogenesis.",
#     metadata = read.table("data/edgeR_example4_GSE60450_Lactation_metadata.tab", header=TRUE)
#     ),
#   tuch = list(
#     dataframe = (function() { 
#       ids <- c("N8", "N33", "N51", "T8", "T33", "T51")
#       filenames <- file.path("data", paste0("Tuch_", ids, ".tab"))
#       tab <- do.call(cbind, lapply(filenames, function(x) read.table(x)[, 2]))
#       
#       colnames(tab) <- ids
#       rownames(tab) <- read.table(filenames[1])[, 1]
# 
#       return(as.data.frame(tab))
#     })(),
#     description = "Tuch et al."
#   )
# )

# example datasets
example.datasets <- list(
  lactate = list(
    dataframe = read.table("data/edgeR_example4_GSE60450_Lactation-GenewiseCounts.tab", header=TRUE, row.names = 1)[, -1],
    description = "Lactogenesis.",
    metadata = read.table("data/edgeR_example4_GSE60450_Lactation_metadata.tab", header=TRUE)
  )
)

ensembl <- useEnsembl(biomart="ensembl")

# ui <- navbarPage(id = "tabs",
#                  theme = "bootstrap.css",
#                  title = "D-fferrential",
#                  tabPanel("Import", tab_importUI("tab_import", example.datasets)),
#                  tabPanel("Visualize", tab_visualizeUI("tab_visualize")),
#                  tabPanel("Differential Expression", tab_deUI("tab_de")),
#                  tabPanel("GO Enrichment", mod_GOUI("mod_GO", ensembl)),
#                  #tabPanel("Additional Plots", tab_plotsUI("tab_plots")),
#                  tabPanel("Export Analysis", tab_exportUI("tab_export")))
                 
ui <- dashboardPage(
  dashboardHeader(title = "D-fferrential"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Import Data", tabName="ImportData", icon = icon("cloud-upload"), selected = TRUE),
      #menuItem("Visualize Data", tabName="VisualizeData"),
      menuItem("Differential Expression", tabName="DifferentialExpression", icon = icon("not-equal")),
      menuItem("GO Enrichment", tabName="GOEnrichment", icon = icon("sitemap")), # aggs: add icon
      menuItem("Export Analysis", tabName = "ExportAnalysis", icon = icon("cloud-download") )) # aggs: add icon
  ),
  dashboardBody(
    tabItems(
      tabItem("ImportData", tab_importUI("tab_import", example.datasets)),
      #tabItem("VisualizeData", tab_visualizeUI("tab_visualize")),
      tabItem("DifferentialExpression", tab_deUI("tab_de")),
      tabItem("GOEnrichment", mod_GOUI("mod_GO", ensembl)),
      tabItem("ExportAnalysis", tab_exportUI("tab_export")))
  )
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
  #sessionData <- callModule(tab_visualizeServer, "tab_visualize", sessionData)
  
  # DE tab
  sessionData <- callModule(tab_deServer, "tab_de", sessionData)
  
  # GO tab
  sessionData <- callModule(mod_GOServer, "mod_GO", sessionData, ensembl)
  
  # plots tab
  #sessionData <- callModule(tab_plotsServer, "tab_plots", sessionData)
  
  # export
  sessionData <- callModule(tab_exportServer, "tab_export", sessionData)

  # observeEvent(sessionData$dataframe(), {
  #   req(sessionData$dataframe())
  #   
  #   showTab("tabs", target = "Filter")
  # })
  
}

shinyApp(ui, server)


