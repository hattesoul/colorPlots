#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# A Forward-Pipe Operator for R
library(magrittr)

# Colored Terminal Output
require(crayon)

# Interpreted String Literals
require(glue)

# Create Interactive Web Graphics via 'plotly.js'
require(plotly)

# HTML Widgets for R
require(htmlwidgets)

# Methods for Single-Cell RNA-Seq Data Analysis
library(scran)

# Single-Cell Analysis Toolkit for Gene Expression Data in R
library(scater)

# Utility functions for data from single-cell experiments
# remotes::install_github("jenzopr/singlecellutils")
library(singlecellutils)

# Functions for timing R scripts, as well as implementations of Stack and List structures
library(tictoc)

# set upload file size limit to 500 GB
options(shiny.maxRequestSize = 500*1024^2)

# load data object
# object <- NULL
object <- readRDS("assets/object_06_compressed.rds")
# saveRDS(object, glue("assets/object_06_compressed.rds"))

# data(diamonds, package = "ggplot2")
nms <- c("Library", "Patient", "Condition", "Total UMI", "Unique features", "(X) Log counts", paste(sort(rowData(object)["symbol"][[1]])))
myOpacity <- NULL

cat("BEGIN\n")

ui <- fluidPage(
  tags$head(tags$script('
        var dimension = [0, 0];
        $(document).on("shiny:connected", function(e) {
          dimension[0] = window.innerWidth;
          dimension[1] = window.innerHeight;
          Shiny.onInputChange("dimension", dimension);
        });
        $(window).resize(function(e) {
           dimension[0] = window.innerWidth;
           dimension[1] = window.innerHeight;
           Shiny.onInputChange("dimension", dimension);
        });
      ')),
  
  # headerPanel("Color plots"),
  sidebarLayout(position = "right",
    sidebarPanel(
      h3("Color plots"),
      selectInput("color", "Select color", choices = nms, selected = "Library"),
      h3("Load cluster file"),
      fileInput("inputFile", "Choose cluster RDS file", placeholder = ".rds file")
      # h3("Save plot file"),
      # downloadButton('downloadPlot', 'Download Plot')
    ),
    mainPanel(
      br(),
      fillPage(plotlyOutput("trendPlot"))
    )
  )
)

server <- function(input, output) {

  output$dimension_display <- renderText({
    paste(input$dimension[1], input$dimension[2], input$dimension[2]/input$dimension[1])
  })

  plotListUMAP2D <- vector("list")
  
  output$trendPlot <- renderPlotly({

    if(is.null(object)) {
    } else {

      # build graph with plotly

      # set colors
      if(match(input$color, nms) > 6) {
        myColor <- assay(object, "umi")[match(input$color, as.list(rowData(object)["symbol"])[[1]]), ]
        myOpacity <- ifelse(assay(object, "umi")[match(input$color, as.list(rowData(object)["symbol"])[[1]]), ] != 0, 0.6, 0.05)
      } else {
        myOpacity <- ifelse(as.list(colData(object)[match("total_umi", names(colData(object)))])[[1]] != 0, 0.6, 0.05)
        if(match(input$color, nms) == 4) {
          myColor <- as.list(colData(object)[match("total_umi", names(colData(object)))])[[1]]
        } else if(match(input$color, nms) == 5) {
          myColor <- as.list(colData(object)[match("total_features_by_umi", names(colData(object)))])[[1]]
        } else if(match(input$color, nms) == 6) {
          myColor <- colSums(assay(object, "logcounts"))
        } else {
          myColor <- as.list(colData(object)[match(input$color, names(colData(object)))])[[1]]
        }
      }
  
      plotNames <- names(object@reducedDims)
      plotListUMAP2D <- vector("list")
  
      p <- plot_ly()
      
      for (i in 1:length(plotNames)) {
        currentParameters = strsplit(plotNames[i], "_")
        myParameterName1 <- NULL
        myParameterName2 <- NULL
        myParameterValue1 <- NULL
        myParameterValue2 <- NULL
  
        myParameterName1 = "minimum distance"
        myParameterName2 = "neighbors"
        if(pmatch("dis", currentParameters[[1]][3], nomatch = 0)) {
          myParameterValue1 <- regmatches(currentParameters[[1]][3], regexpr("\\d+\\.?(\\d+)?", currentParameters[[1]][3], perl = TRUE))
        }
        if(pmatch("nei", currentParameters[[1]][4], nomatch = 0)) {
          myParameterValue2 <- regmatches(currentParameters[[1]][4], regexpr("\\d+", currentParameters[[1]][4], perl = TRUE))
        }
        
        xlabel <- list(title = glue(myParameterName2, ": ", myParameterValue2))
        ylabel <- list(title = glue(myParameterName1, ": ", myParameterValue1))
        
        # p <- plot_ly(alpha = 0.6, x = reducedDim(object, plotNames[i])[, 1], y = reducedDim(object, plotNames[i])[, 2], type = "scattergl", mode = "markers", color = myColor, legendgroup = myColor, text = myColor, marker = list(opacity = myOpacity, colorscale = "Viridis", reversescale = TRUE, showscale = FALSE), showlegend = ifelse(i == 1, TRUE, FALSE), height = round((as.numeric(input$dimension[2]) * 0.95))) %>%
        p <- plot_ly(alpha = 0.6, x = reducedDim(object, plotNames[i])[, 1], y = reducedDim(object, plotNames[i])[, 2], type = "scattergl", mode = "markers", color = myColor, legendgroup = myColor, text = myColor, colorbar = NULL, marker = list(opacity = myOpacity, reversescale = TRUE, showscale = ifelse(i == 1, FALSE, FALSE)), showlegend = ifelse(i == 1, TRUE, FALSE), height = round((as.numeric(input$dimension[2]) * 0.95))) %>%
          layout(title = "", xaxis = xlabel, yaxis = ylabel)
  
        plotListUMAP2D[[i]] <- p
      }
      # p
      plotsUMAP2D <- subplot(plotListUMAP2D, nrows = 3, titleX = TRUE, titleY = TRUE, shareX = TRUE, shareY = TRUE) %>%
        layout(title = "UMAP")
    }
  })

  # download plot file
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(glue("{inputFile}_{input$color}.html"))
    },
    content = function(file) {
      # saveWidgetFix(as_widget(plotsUMAP2D), savefile)
      myWidget <- htmlwidgets::saveWidget(plotsUMAP2D, savefile)
      # write.csv(datasetInput(), file, row.names = FALSE)
      myWidget
    }
  )
}

shinyApp(ui = ui, server = server)
