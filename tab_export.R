
#' UI function for statistics module
tab_exportUI <- function(id) {
  ns <- NS(id)
  
  # summary.ui <- tagList(
  #   h3("Summary"),
  #   htmlOutput(ns("text_summary"))
  # )
  
  tagList(
    useShinyjs(),
    box(width=12,
      h4("Differential expression results"),
      uiOutput(ns("ui_de_tables"))),
    box(width=12,
      h4("GO enrichment results"),
      uiOutput(ns("ui_go_tables"))),
    box(width=12,
      verbatimTextOutput(ns("txt_session")))
    #h4("Analysis report"),
    #inputPanel(radioButtons(ns("radio_format"), label = "Output format", choices = c("HTML", "Pdf", "Doc"))),
    
    #summary.ui,
    #uiOutput(ns("uiSections")),
    
    #downloadButton(ns("report"), "Generate report")
  )
}

#' Server function for statistics module
#' 
#' @param cmatrix Counts matrix.
#' 
#' @return A dataframe as a reactive value.
tab_exportServer <- function(input, output, session, sessionData) {
  
  # output$text_summary <- renderUI({
  #   print(sessionData)
  #   
  #   fparams <- sessionData$filter_params()
  #   
  #   req(fparams)
  #   
  #   ftext <- tagList(
  #     h4("Filtering parameters"),
  #     tags$ul(
  #       tags$li("Minimum UMI per cell: ", fparams$min_umi),
  #       tags$li("Maximum UMI per cell: ", fparams$max_umi),
  #       tags$li("Minimum genes per cell: ", fparams$min_genes),
  #       tags$li("Maximum genes per cell: ", fparams$max_genes))
  #   )
  #   
  #   return(ftext)
  # })
  
  # disable("report")
  
  # observe({
  #   req(sessionData$all_markers())
  # 
  #   enable("report")
  # })
  
  create.download.button <- function(id, label, tab) {
    button <- downloadButton(session$ns(id), label = label)
    
    # create the downloadHandler for the button
    output[[ id ]] <- downloadHandler(
      filename = function() {
        paste0("d-fferential-", Sys.Date(), "-", id, ".csv")
      },
      content = function(file) {
        write.csv(tab, file)
      }
    )
    
    return(button)
  }
  
  output$ui_de_tables <- renderUI({
    data <- sessionData$de_results$data
    
    header <- tags$tr(tags$th("ID"), tags$th("Method"), tags$th("FDR cutoff"), tags$th("logFC cutoff"), tags$th("Results"))
    
    rows <- lapply(data, function(x) {
      button <- create.download.button(paste0("download_de_csv_", x$id), "csv", x$tab)
      
      # make the row
      tags$tr(tags$td(x$id), tags$td(x$method), tags$td(x$fdr), tags$td(x$lfc), tags$td(button))
    })
    
    htmltab <- tags$table(header, rows, class="table")
    
    return(htmltab)
  })

  output$ui_go_tables <- renderUI({
    data <- sessionData$go_results
    
    header <- tags$tr(tags$th("ID"), tags$th("Geneset source"), tags$th("Geneset"), tags$th("Annotation"), tags$th("Gene IDs"),
                      tags$th("BP"), tags$th("CC"), tags$th("MF"))
    rows <- lapply(data, function(x) {
      button.bp <- create.download.button(paste0("download_go_bp_csv_", x$id), "csv", x$bp.tab)
      button.cc <- create.download.button(paste0("download_go_cc_csv_", x$id), "csv", x$cc.tab)
      button.mf <- create.download.button(paste0("download_go_mf_csv_", x$id), "csv", x$mf.tab)
      
      # make the row
      tags$tr(tags$td(x$id), tags$td(x$geneset.source), tags$td(x$geneset.type), tags$td(x$annotation.name), tags$td(x$attribute.name),
              tags$td(button.bp), tags$td(button.cc), tags$td(button.mf))
    })
    
    htmltab <- tags$table(header, rows, class="table")
    
    return(htmltab)
  })
  
  
  output$report <- downloadHandler(
    filename = function() {
      paste0("d-fferential-report-", Sys.Date(), ".html")
    },
    content = function(file) {
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      plots <- lapply(sessionData$plots, function(x) x$plot.fun())
      
      # Set up parameters to pass to Rmd document
      params <- list(dataframe = sessionData$dataframe(),
                     metadata = sessionData$metadata(),
                     plots = plots,
                     de_results = sessionData$de_results)
      
      withProgress(message = "Generating report...", {
        rmarkdown::render(tempReport,
                          output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv()))
      })
    }
  )
  
  output$txt_session <- renderPrint({
    sessionInfo()
  })
  
  return(sessionData)
  
}

