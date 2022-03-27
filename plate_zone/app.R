library(rhandsontable)
library(rjags)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(tidybayes)
library(shiny)
theme_set(theme_cowplot(font_size = 18))

ui <- shinyUI(fluidPage(
  title = "Analysis of zones of inhibition",
  titlePanel("Analysis of zones of inhibition"),
  fluidPage(
    column(
      width = 4,
      h3("Experimental data"),
      helpText("Use this table to enter your own experimental data. You can edit each cell, add and remove rows with right click. Enter the diameter of the zone of inhibition in millimeter, a unique id per disk treatment, a unique number for each plate you are recording, and a unique id for each student. Once you are done, click Update plots."),
      rHandsontableOutput("hot"),
      br(),
      submitButton(text = "Update plots", icon("refresh"), width = "100%")
    ),
    column(
      width = 8,
      tabsetPanel(type = "tabs",
                  tabPanel("Plot",
                           helpText("This plot shows your data unprocessed. Students are color coded and plates are shape coded. Variations in the zone diameters between plates and students may make it difficult to determine if the different treatments have different effects. A statistical analysis is required to determine how confident we are that the treatments have different effects on growth"),
                           plotOutput("plotRawData",height = "400px")
                  ),
                  tabPanel("Analysis",
                           helpText("The first plot shows the credible probability distributions of the zones of inhibitions for each treatment given the data and the statistical model. Because of experimental errors and biological variability, the credible intervals can be wide. Increased experimental precision and replication can help narrow the distributions by providing more information to the analysis."),
                           plotOutput("plotPost",height = "300px"),
                           helpText("The second plot shows the credible probability distributions of the difference between each pair of treatments. This information can be used to determine how significant these differences are, but also, to evaluate the magnitudes of the differences."),
                           plotOutput("plotDiffDisk",height = "300px"),
                           helpText("The last two plots can help evaluate if the variability in the data can be attributed to other experimental factors, such as, differences in plate preparation or student technique."),
                           plotOutput("plotDiffPlateStudent",height = "300px")
                  ),
                  tabPanel("Model",
                           tags$div(class="header", checked=NA,
                                    br(),
                                    tags$h5("The data was analyzed using a Bayesian approach."),
                                    tags$p("The diameters of the zones of inhibition are modeled as linear mixture of non-negative treatment effects with random effects from plates and students."),
                                    tags$h5("Model"),
                                    tags$p("Zone ~ (1|Disk) + (1|Plate) + (1|Student)"),
                                    tags$h5("Priors"),
                                    tags$p("Disk ~ Exponential(TauD)"),
                                    tags$p("Plate ~ Normal(0,SigP)"),
                                    tags$p("Student ~ Normal(0,SigS)")
                           )
                  )
      )
    )
  )
))

server <- shinyServer(
  function(input, output) {

    DF <-
      data.frame(
        Zone = c(0, 10, 15, 7, 1, 10.5, 11.5, 11, 3, 15, 21, 14),
        Disk = c("A", "B", "C", "D", "A", "B", "C", "D", "A", "B", "C", "D"),
        Plate = c("1", "1", "1", "1", "2", "2", "2", "2", "3", "3", "3", "3"),
        Student = c(
          "JJ",
          "JJ",
          "JJ",
          "JJ",
          "JJ",
          "JJ",
          "JJ",
          "JJ",
          "CC",
          "CC",
          "CC",
          "CC"
        ),
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
      DF$Disk = as.factor(DF$Disk)
      DF$Plate = as.factor(DF$Plate)
      DF$Student = as.factor(DF$Student)
      p <-
        ggplot(DF) + geom_point(aes(
          y = fct_rev(Disk),
          x = Zone,
          shape = Plate,
          color = Student
        ), size = 6) + expand_limits(y = 0) + xlab('Zone diameter (mm)') + ylab('Disk')
      print(p)
    })

    samples <- reactive({
      withProgress(message = 'Calculation in progress...', value = 0, {
        DF <- values[["DF"]]
        DF$Disk = as.factor(DF$Disk)
        DF$Plate = as.factor(DF$Plate)
        DF$Student = as.factor(DF$Student)

        data = compose_data(DF)

        print(data)

        params_inits = list(
          "sigma" = 10,
          "sigma_d" = 10,
          "sigma_p" = 10,
          "sigma_s" = 10
        )

        model <-
          jags.model(
            "plate_zone.bug",
            data = data,
            inits = params_inits,
            n.chains = 4
          )

        update(model, n.iter = 10000)

        samples <-
          coda.samples(
            model,
            variable.names = c("d", "p", "s"),
            n.iter = 2500,
            thin = 10
          )

        samples = recover_types(samples, DF)

      })
    })

    output$plotPost = renderPlot({
      samples = samples()
      p_post_disk = samples %>%
        spread_draws(d[Disk]) %>%
        ggplot(aes(y = fct_rev(Disk), x = d)) + geom_eyeh() + xlab('Zone (mm)') + ylab('Disk')
      print(p_post_disk)
    })

    output$plotDiffDisk = renderPlot({
      samples = samples()
      p_diff_disk = samples %>%
        spread_draws(d[Disk]) %>%
        compare_levels(d, by = Disk) %>%
        ggplot(aes(y = fct_rev(Disk), x = d)) + geom_eyeh() + geom_vline(xintercept = 0) + xlab('Zone difference (mm)') + ylab('Disk')
      print(p_diff_disk)
    })
    output$plotDiffPlateStudent = renderPlot({
      samples = samples()
      p_diff_plate = samples %>%
        spread_draws(p[Plate]) %>%
        compare_levels(p, by = Plate) %>%
        ggplot(aes(y = fct_rev(Plate), x = p)) + geom_eyeh() + geom_vline(xintercept = 0) + xlab('Zone difference (mm)') + ylab('Plate')

      p_diff_student = samples %>%
        spread_draws(s[Student]) %>%
        compare_levels(s, by = Student) %>%
        ggplot(aes(y = fct_rev(Student), x = s)) + geom_eyeh() + geom_vline(xintercept = 0)  + xlab('Zone difference (mm)') + ylab('Student')

      p = plot_grid(p_diff_plate, p_diff_student, ncol = 2)
      print(p)
    })
  })
## run app
shinyApp(ui = ui, server = server)
