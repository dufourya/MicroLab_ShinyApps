library(rhandsontable)
library(rjags)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(tidybayes)
library(shiny)
theme_set(theme_cowplot(font_size = 18))

ui <- shinyUI(fluidPage(
    title = "Comparing groups",
    titlePanel("Comparing groups"),
    fluidPage(
        column(
            width = 4,
            h3("Experimental data"),
            helpText("Use this table to enter your own experimental data. You can edit each cell, add and remove rows with right click. Label the groups you would like to compare and the nuisance factors that may confound your analysis."),
            rHandsontableOutput("hot"),
            br(),
            submitButton(text = "Update plots", icon("refresh"), width = "100%")
        ),
        column(
            width = 8,
            tabsetPanel(type = "tabs",
                        tabPanel("Plot",
                                 helpText("This plot shows your data unprocessed. The nuisance factors are indicated with shapes and colors. The contributions of the nuisance factors to the data makes it difficult to determine if the difference between the groups is significant. A statistical analysis is required to determine how confident we are that the groups have different means."),
                                 plotOutput("plotRawData",height = "400px")
                        ),
                        tabPanel("Analysis",
                                 helpText("The first plot shows the credible probability distributions of the means for each group given the data and the model. Because of random effects from measurements and nuisance factors, the credible intervals can be wide. Increased precision and replication can help narrow the distributions by providing more information for the analysis."),
                                 plotOutput("plotPost",height = "300px"),
                                 helpText("The second plot shows the credible probability distributions of the difference between the group means without the nuisance factors. This information can be used to determine how significant the difference is, but also, to evaluate the magnitude of the difference."),
                                 plotOutput("plotDiffGroup",height = "300px"),
                                 helpText("The last two plots can help evaluate if the variability in the data can be attributed to other nuisance factors."),
                                 plotOutput("plotDiffFactor1Factor2",height = "300px")
                        ),
                        tabPanel("Model",
                                 tags$div(class="header", checked=NA,
                                          br(),
                                          tags$h5("The data was analyzed using a Bayesian approach."),
                                          tags$p("The data are modeled the sum of group effects with random effects from Factor1 and Factor2."),
                                          tags$h5("Model"),
                                          tags$p("Data ~ (1|Group) + (1|Factor1) + (1|Factor2)"),
                                          tags$h5("Priors"),
                                          tags$p("Group ~ Normal(meanData,SigD)"),
                                          tags$p("Factor1 ~ Normal(0,SigP)"),
                                          tags$p("Factor2 ~ Normal(0,SigS)")
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
                Data = c(rnorm(10, 10, 1), rnorm(10, 13, 2)),
                Group = c(rep("A",10), rep("B",10)),
                Factor1 = sample(c(rep("JJ",12), rep("CC",8))),
                Factor2 = sample(c(rep("1",7), rep("2",7), rep("3",6))),
                stringsAsFactors = FALSE
            )
        DF$Data[DF$Factor1 == "JJ"] = DF$Data[DF$Factor1 == "JJ"] + rnorm(12,1,0.1)

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
            DF$Group = as.factor(DF$Group)
            DF$Factor1 = as.factor(DF$Factor1)
            DF$Factor2 = as.factor(DF$Factor2)
            p <-
                ggplot(DF) + geom_point(aes(
                    y = fct_rev(Group),
                    x = Data,
                    shape = Factor1,
                    color = Factor2
                ), size = 6) + expand_limits(y = 0) + xlab('Data') + ylab('Group')
            print(p)
        })

        samples <- reactive({
            withProgress(message = 'Calculation in progress...', value = 0, {
                DF <- values[["DF"]]
                DF$Group = as.factor(DF$Group)
                DF$Factor1 = as.factor(DF$Factor1)
                DF$Factor2 = as.factor(DF$Factor2)

                data = compose_data(DF)
                data$meanData = mean(data$Data)

                #print(data)

                params_inits = list(
                    "sigma" = 10,
                    "sigma_d" = 10,
                    "sigma_f1" = 10,
                    "sigma_f2" = 10
                )

                model <-
                    jags.model(
                        "comparing_groups.bug",
                        data = data,
                        inits = params_inits,
                        n.chains = 4
                    )

                update(model, n.iter = 10000)

                samples <-
                    coda.samples(
                        model,
                        variable.names = c("d", "f1", "f2"),
                        n.iter = 5000,
                        thin = 50
                    )

                samples = recover_types(samples, DF)

            })
        })

        output$plotPost = renderPlot({
            samples = samples()
            p_post_Group = samples %>%
                spread_draws(d[Group]) %>%
                ggplot(aes(y = fct_rev(Group), x = d)) + geom_eyeh(fill = 'firebrick') + xlab('Data') + ylab('Group') + geom_vline(xintercept = 0, linetype = 'dashed')
            print(p_post_Group)
        })

        output$plotDiffGroup = renderPlot({
            samples = samples()
            p_diff_Group = samples %>%
                spread_draws(d[Group]) %>%
                compare_levels(d, by = Group) %>%
                ggplot(aes(y = fct_rev(Group), x = d)) + geom_eyeh(fill = 'firebrick') + geom_vline(xintercept = 0, linetype = 'dashed') + ylab('Group')
            print(p_diff_Group)
        })
        output$plotDiffFactor1Factor2 = renderPlot({
            samples = samples()
            p_diff_Factor1 = samples %>%
                spread_draws(f1[Factor1]) %>%
                compare_levels(f1, by = Factor1) %>%
                ggplot(aes(y = fct_rev(Factor1), x = f1)) + geom_eyeh(fill = 'firebrick') + geom_vline(xintercept = 0, linetype = 'dashed') + xlab('Difference') + ylab('Factor1')

            p_diff_Factor2 = samples %>%
                spread_draws(f2[Factor2]) %>%
                compare_levels(f2, by = Factor2) %>%
                ggplot(aes(y = fct_rev(Factor2), x = f2)) + geom_eyeh(fill = 'firebrick') + geom_vline(xintercept = 0, linetype = 'dashed')  + xlab('Difference') + ylab('Factor2')

            p = plot_grid(p_diff_Factor1, p_diff_Factor2, ncol = 2)
            print(p)
        })
    })
## run app
shinyApp(ui = ui, server = server)
