library(shiny)
library(ggplot2)
library(cowplot)
library(rjags)
#library(tidyverse)
#library(tidybayes)
theme_set(theme_cowplot())

ui <- fluidPage(
  titlePanel("Logistic growth"),
  sidebarLayout(
    sidebarPanel(
      sliderInput('time_points', 'Sampling interval (min)', min = 5, max = 30, value = 15, step = 1, round = TRUE),
      sliderInput('initial_OD', 'Initial culture OD', min = 0, max = 1, value = 0.1, step = 0.05, round = FALSE),
      sliderInput('growth_rate', 'Generation time (min)', min = 15, max = 75, value = 40, step = 1, round = FALSE),
      sliderInput('carrying_capacity', 'Carrying capacity', min = 0, max = 1, value = 0.5, step = 0.05, round = FALSE),
      sliderInput('proportion_dead_cells', 'Proportion of dead cells', min = 0, max = 1, value = 0.1, step = 0.05,round = FALSE),
      sliderInput('experimental_precision', 'Experimental precision', min = 5, max = 10, value = 8, step = 0.5, round = FALSE),
      #sliderInput('biological_variability', 'Biological variability', min = 0, max = 0.2, value = 0.1, step = 0.01,round = FALSE),
      hr(),
      radioButtons('y_scale', 'Y-axis scaling', choices = list('linear' = 0,'logarithmic' = 1), selected = 0),
      #radioButtons('analysis', 'Data analysis', choices = list('Data averaging' = 0,'Independent fitting' = 1, "Hierarchical modeling" = 2), selected = 0),
      actionButton("fit_button", "Fit data")
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Growth curve",
                           h3("Experimental samples of cell culture optical density"),
                           plotOutput('plot_growth_curve', height = "600px")),
                  tabPanel("Parameter estimates",
                           h3("Posterior probability densities of model parameters"),
                           plotOutput("plot_fitted_parameters", height = "600px")),
                  tabPanel("Parameters scatter",
                           h3("Scatter plots of model parameters"),
                           plotOutput("plot_scatter_parameters", height = "600px"))
      )
    )
  )
)


p = ggplot()

server <- function(input, output, session) {
  session$onSessionEnded(stopApp)
  time_points_fine = seq(0, 500, 1)

  time_points = reactive({
    time_points_fine[seq(1, 501, input$time_points)]
  })

  OD = reactive({
    input$initial_OD * input$proportion_dead_cells + input$carrying_capacity * 2^(time_points_fine / input$growth_rate) / (
      2^(time_points_fine/input$growth_rate) - 1 + input$carrying_capacity / (input$initial_OD * (1 - input$proportion_dead_cells)))
  })

  OD_noisy = reactive({
    OD = OD()
    OD[seq(1, 501, input$time_points)] + rnorm(length(time_points()), 0, sqrt(1 / exp(input$experimental_precision)))
  })

  output$plot_growth_curve <- renderPlot({

    p = p + geom_point(aes(x = time_points(), y = OD_noisy()), size = 2)
    p = p + xlim(0, 500) + xlab("Time (min)") + ylab("Optical density")

    if (input$y_scale == 1) {
      p = p + scale_y_continuous(trans = 'log10', limits = c(0.03, max(1,OD_noisy())))
    } else {
      p = p + scale_y_continuous(trans = 'identity', limits = c(min(0,OD_noisy()), max(1,OD_noisy())))
    }

    if (input$fit_button > 0) {
      OD_fit_all = data.frame(OD = double(), Time = double(), id = integer())
      samples_data = samples()
      for (i in seq(1,1000,10)) {
        OD_fit = data.frame(OD = samples_data$initial_OD[i] * samples_data$proportion_dead_cells[i] + samples_data$carrying_capacity[i] * 2^(time_points_fine/samples_data$growth_rate[i]) / (2^(time_points_fine/samples_data$growth_rate[i]) - 1 + samples_data$carrying_capacity[i] / (samples_data$initial_OD[i] * (1 - samples_data$proportion_dead_cells[i]))), Time = time_points_fine, id = rep(i,length(time_points_fine)))
        OD_fit_all = rbind(OD_fit_all, OD_fit)
      }

      p = p + geom_line(data = OD_fit_all, aes(x = Time, y = OD, group = id), alpha = 0.05)

    }
    
    p = p + geom_line(aes(x = time_points_fine, y = OD()), color = 'firebrick', size = 1)
    
    print(p)
  })
 
  samples <- eventReactive(input$fit_button, {
    
    withProgress(message = 'Fitting parameters', {
      
      growth_curve = list(
        "time_points" = time_points(),
        "OD" = OD_noisy(),
        "N" = length(OD_noisy())
      )
      params_inits = list(
        "experimental_precision" = input$experimental_precision,
        "initial_OD" = input$initial_OD,
        "proportion_dead_cells" = input$proportion_dead_cells,
        "growth_rate" = input$growth_rate,
        "carrying_capacity" = input$carrying_capacity
      )
      
      model <-
        jags.model(
          "growth_model.bug",
          data = growth_curve,
          inits = params_inits,
          n.chains = 4
        )
      
      update(model, n.iter = 10000)
      samples <-
        coda.samples(
          model,
          variable.names = c(
            "experimental_precision",
            "initial_OD",
            "proportion_dead_cells",
            "growth_rate",
            "carrying_capacity"
          ),
          n.iter = 2500,
          thin = 10
        )
      
      samples = do.call(rbind.data.frame, samples)
      
    })
    

  })

  output$plot_fitted_parameters <- renderPlot({

    samples_data = samples()

    p_growth_rate = ggplot(samples_data, aes(x = growth_rate)) + geom_density(fill = "skyblue1") + xlab('Generation time (min)')  + geom_vline(xintercept = input$growth_rate, color = "firebrick") + xlim(15,75)
    p_initial_OD = ggplot(samples_data, aes(x = initial_OD)) + geom_density(fill = "skyblue1") + xlab('Initial OD') + geom_vline(xintercept = input$initial_OD, color = "firebrick") + xlim(0, 1)
    p_carrying_capacity = ggplot(samples_data, aes(x = carrying_capacity)) + geom_density(fill = "skyblue1") + xlab('Carrying capacity') + geom_vline(xintercept = input$carrying_capacity, color = "firebrick") + xlim(0, 1)
    p_proportion_dead_cells = ggplot(samples_data, aes(x = proportion_dead_cells)) + geom_density(fill = "skyblue1") + xlab('Proportion of dead cells') + geom_vline(xintercept = input$proportion_dead_cells, color = "firebrick") + xlim(0, 1)
    #p_experimental_precision = ggplot(samples_data, aes(x = experimental_precision)) + geom_density(fill = "skyblue1") + theme(aspect.ratio = 1) + xlab('Exp. precision') + geom_vline(xintercept = input$experimental_precision, color = "firebrick")#+ xlim(5, 10)
    pout = plot_grid(
      p_initial_OD,
      p_growth_rate,
      p_carrying_capacity,
      p_proportion_dead_cells,
      #p_experimental_precision,
      nrow = 2,
      ncol = 2
    )

    print(pout)

  })

  output$plot_scatter_parameters <- renderPlot({

    samples_data = samples()

    p_gr_init = ggplot(samples_data, aes(x = growth_rate, y = initial_OD)) + geom_point(color = "skyblue1", alpha = 0.1) +
      xlab('Generation time')  + ylab("Initial OD") + geom_point(aes(x = input$growth_rate, y = input$initial_OD), color = "firebrick") + theme(aspect.ratio = 1)

    p_gr_carr = ggplot(samples_data, aes(x = growth_rate, y = carrying_capacity)) + geom_point(color = "skyblue1", alpha = 0.1) +
      xlab('Generation time')  + ylab("Carrying cap.") + geom_point(aes(x = input$growth_rate, y = input$carrying_capacity), color = "firebrick") + theme(aspect.ratio = 1)

    p_gr_prop = ggplot(samples_data, aes(x = growth_rate, y = proportion_dead_cells)) + geom_point(color = "skyblue1", alpha = 0.1) +
      xlab('Generation time')  + ylab("Prop. dead cells") + geom_point(aes(x = input$growth_rate, y = input$proportion_dead_cells), color = "firebrick") + theme(aspect.ratio = 1)

    p_init_carr = ggplot(samples_data, aes(x = initial_OD, y = carrying_capacity)) + geom_point(color = "skyblue1", alpha = 0.1) +
      xlab('Initial OD')  + ylab("Carrying cap.") + geom_point(aes(x = input$initial_OD, y = input$carrying_capacity), color = "firebrick") + theme(aspect.ratio = 1)

    p_init_prop = ggplot(samples_data, aes(x = initial_OD, y = proportion_dead_cells)) + geom_point(color = "skyblue1", alpha = 0.1) +
      xlab('Initial OD')  + ylab("Prop. of dead cells") + geom_point(aes(x = input$initial_OD, y = input$proportion_dead_cells), color = "firebrick") + theme(aspect.ratio = 1)

    p_carr_prop = ggplot(samples_data, aes(x = carrying_capacity, y = proportion_dead_cells)) + geom_point(color = "skyblue1", alpha = 0.1) +
      xlab('Carrying capacity')  + ylab("Prop. of dead cells") + geom_point(aes(x = input$carrying_capacity, y = input$proportion_dead_cells), color = "firebrick") + theme(aspect.ratio = 1)

    pout = plot_grid(
      p_gr_init,
      p_gr_carr,
      p_gr_prop,
      p_init_carr,
      p_init_prop,
      p_carr_prop,
      nrow = 3,
      ncol = 2
    )

    print(pout)

  })
}

shinyApp(ui = ui, server = server)
