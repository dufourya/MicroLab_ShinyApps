library(rhandsontable)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(shiny)
theme_set(theme_cowplot(font_size = 18))

ui <- shinyUI(fluidPage(
  title = "Analysis of bacterial growth",
  titlePanel("Analysis of bacterial growth"),
  fluidPage(
    column(
      width = 4,
      h3("Experimental data"),
      helpText("Use this table to enter your own experimental data. You can edit each cell, add and remove rows with right click. Enter the simulation time steps and the corresponding number of cells for each simulation trial. Once you are done, click Update plots."),
      rHandsontableOutput("hot"),
      br(),
      submitButton(text = "Update plots", icon("refresh"), width = "100%")
    ),
    column(
      width = 8,
      tabsetPanel(type = "tabs",
                  tabPanel("Plot",
                           helpText("These plots show the data unprocessed. The top plot shows the cell numbers on a linear scale while the bottom plot uses a logarithmic scale."),
                           plotOutput("plotRawData",height = "800px")
                  ),
                  tabPanel("Analysis",
                           helpText("The first plot shows the growth rate of each population as a function of time. The second plot shows growth rate as a function of cell number on a log scale."),
                           plotOutput("plotAnalysis",height = "800px")
                  )
      )
    )
  )
))

server <- shinyServer(

  function(input, output, session) {

    session$onSessionEnded(function() {
      stopApp()
    })

    DF <-
      data.frame(
        Time = as.integer(c(0, 100, 200, 300, 400, 500, 600, 700, 800, 0, 100, 200, 300, 400, 500, 600, 700, 800, 0, 100, 200, 300, 400, 500, 600, 700, 800)),
        Cells = as.integer(c(1, 7, 35, 123, 218, 364, 531, 716, 892, 1, 4, 14, 57, 207, 484, 854, 1311, 1836, 1, 2, 4, 8, 16, 32, 60, 110, 200)),
        Trial = as.integer(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3)),
        stringsAsFactors = FALSE
      )

    ## Handsontable
    values <- reactiveValues()

    observe({
      if (!is.null(input$hot)) {
        DF = hot_to_r(input$hot)
      } else {
        if (is.null(values[["DF"]]))
          DF <- DF
        else
          DF <- values[["DF"]]
      }
      values[["DF"]] <- DF
    })

    output$hot <- renderRHandsontable({
      DF <- values[["DF"]]
      if (!is.null(DF))
        rhandsontable(
          DF,
          useTypes = as.logical(TRUE),
          stretchH = "all",
          readOnly = FALSE
        )
    })

    output$plotRawData <- renderPlot({
      DF <- values[["DF"]]
      DF$Trial = as.factor(DF$Trial)

      p_lin <-
        ggplot(DF) + geom_point(aes(
          x = Time,
          y = Cells,
          color = Trial
        ), size = 5) + expand_limits(y = 0) + xlab('Time step') + ylab('Number of cells') + theme(legend.position="top")
      p_log <-
        ggplot(DF) + geom_point(aes(
          x = Time,
          y = Cells,
          color = Trial
        ), size = 5) + xlab('Time step') + ylab('Number of cells') + scale_y_continuous(trans ='log2') + theme(legend.position="none")

      p = plot_grid(p_lin, p_log, ncol = 1)
      print(p)
    })

    output$plotAnalysis = renderPlot({
      DF <- values[["DF"]]
      DF$Trial = as.factor(DF$Trial)

      DF = DF %>% group_by(Trial)
      DF = DF %>% mutate(DiffTime = c(diff(Time), NA), DiffCells = c(diff(log2(Cells)),NA)) %>% drop_na()

      p_time <-
        ggplot(DF) + geom_point(aes(
          x = Time + 50,
          y = DiffCells / DiffTime * 50,
          color = Trial
        ), size = 5) + expand_limits(y = 0) + xlab('Time step') + ylab('Growth Rate') + theme(legend.position="top")

      p_cells <-
        ggplot(DF) + geom_point(aes(
          x = Cells + DiffCells/2,
          y = DiffCells / DiffTime * 50,
          color = Trial
        ), size = 5) + expand_limits(y = 0) + xlab('Cell number') + ylab('Growth Rate') + theme(legend.position="none") + scale_x_continuous(trans ='log2')

      p = plot_grid(p_time, p_cells, ncol = 1)
      print(p)
    })

  })
## run app
shinyApp(ui = ui, server = server)
