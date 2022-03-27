library(rhandsontable)
library(rjags)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(modelr)
library(tidybayes)
library(shiny)
theme_set(theme_cowplot(font_size = 18))

ui <- shinyUI(fluidPage(
    title = "Linear regression",
    titlePanel("Linear regression"),
    fluidPage(
        column(
            width = 4,
            h3("Experimental Variable"),
            helpText("Use this table to enter your own experimental Variable. You can edit each cell, add and remove rows with right click. Label the Responses you would like to compare and the nuisance factors that may confound your analysis."),
            rHandsontableOutput("hot"),
            br(),
            submitButton(text = "Update plots", icon("refresh"), width = "100%")
        ),
        column(
            width = 8,
            tabsetPanel(type = "tabs",
                        tabPanel("Plot",
                                 helpText("This plot shows your Variable unprocessed. The nuisance factors are indicated with shapes and colors. The contributions of the nuisance factors to the Variable makes it difficult to determine if the difference between the Responses is significant. A statistical analysis is required to determine how confident we are that the Responses have different means."),
                                 plotOutput("plotRawVariable",height = "400px")
                        ),
                        tabPanel("Analysis",
                                 helpText("The first plot shows the credible probability distributions of the means for each Response given the Variable and the model. Because of random effects from measurements and nuisance factors, the credible intervals can be wide. Increased precision and replication can help narrow the distributions by providing more information for the analysis."),
                                 plotOutput("plotRegression",height = "400px"),
                                 helpText("The second plot shows the credible probability distributions of the difference between the Response means without the nuisance factors. This information can be used to determine how significant the difference is, but also, to evaluate the magnitude of the difference."),
                                 plotOutput("plotPost",height = "300px"),
                                 helpText("The last two plots can help evaluate if the variability in the Variable can be attributed to other nuisance factors."),
                                 plotOutput("plotDiffResponse",height = "300px")
                        ),
                        tabPanel("Model",
                                 tags$div(class="header", checked=NA,
                                          br(),
                                          tags$h5("The Variable was analyzed using a Bayesian approach."),
                                          tags$p("The Variable are modeled the sum of Response effects with random effects from Group and Factor2."),
                                          tags$h5("Model"),
                                          tags$p("Variable ~ (1|Response) + (1|Group) + (1|Factor2)"),
                                          tags$h5("Priors"),
                                          tags$p("Response ~ Normal(meanVariable,SigD)"),
                                          tags$p("Group ~ Normal(0,SigP)"),
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
                Variable = round(c(seq(1,20, 2), seq(3,17,1.5)),digits = 1),
                Response = round(c(seq(1,20, 2) * 1.5, seq(3,17,1.5) * 2.25) + c(rnorm(10,15,2.25), rnorm(10,10,2.25)), digits = 1),
                Group = c(rep("JJ",10), rep("CC",10)),
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

        output$plotRawVariable <- renderPlot({
            DF <- values[["DF"]]
            DF$Group = as.factor(DF$Group)
            p <-
                ggplot(DF) + geom_point(aes(
                    y = Response,
                    x = Variable,
                    color = Group
                ), size = 3) + xlab('Variable') + ylab('Response') + expand_limits(y = c(10,50), x = c(0,21))
            print(p)
        })

        samples <- reactive({
            withProgress(message = 'Calculation in progress...', value = 0, {
                DF <- values[["DF"]]
                DF$Group = as.factor(DF$Group)

                data = compose_data(DF)

                #print(Variable)

                params_inits = list(
                    "sigma" = 10,
                    "sigma_a" = 10,
                    "sigma_b" = 10
                )

                model <-
                    jags.model(
                        "linear_regression.bug",
                        data = data,
                        inits = params_inits,
                        n.chains = 4
                    )

                update(model, n.iter = 10000)

                samples <-
                    coda.samples(
                        model,
                        variable.names = c("a", "b"),
                        n.iter = 2500,
                        thin = 10
                    )

                samples = recover_types(samples, DF)

            })
        })

        output$plotPost = renderPlot({
            samples = samples()
            p_post_Slope = samples %>%
                spread_draws(a[Group]) %>%
                ggplot(aes(y = Group, x = a)) + geom_eyeh(fill = 'firebrick') + xlab('Slope') + ylab('Group') + geom_vline(xintercept = 0, linetype = 'dashed')
            p_post_Intercept = samples %>%
                spread_draws(b[Group]) %>%
                ggplot(aes(y = Group, x = b)) + geom_eyeh(fill = 'firebrick') + xlab('Intercept') + ylab('Group') + geom_vline(xintercept = 0, linetype = 'dashed')

            p = plot_grid(p_post_Slope, p_post_Intercept, ncol = 2)
            print(p)
        })

        output$plotDiffResponse = renderPlot({
            samples = samples()
            p_diff_Slope = samples %>%
                spread_draws(a[Group]) %>%
                compare_levels(a, by = Group) %>%
                ggplot(aes(y = Group, x = a)) + geom_eyeh(fill = 'firebrick') + geom_vline(xintercept = 0, linetype = 'dashed') + xlab('Difference') + ylab('Slope')
            p_diff_Intercept = samples %>%
                spread_draws(b[Group]) %>%
                compare_levels(b, by = Group) %>%
                ggplot(aes(y = Group, x = b)) + geom_eyeh(fill = 'firebrick') + geom_vline(xintercept = 0, linetype = 'dashed') + xlab('Difference') + ylab('Intercept')

            p = plot_grid(p_diff_Slope, p_diff_Intercept, ncol = 2)
            print(p)
        })

        output$plotRegression = renderPlot({

            DF <- values[["DF"]]
            DF$Group = as.factor(DF$Group)

            samples <- samples()
            ablines = samples %>%
                spread_draws(a[Group], b[Group])

            p_Regression = ggplot(DF) + geom_point(aes(x = Variable, y = Response, color = Group), size = 3) +
                expand_limits(y = c(10,50), x = c(0,21)) + xlab('Variable') + ylab('Response')
            p_Regression = p_Regression + geom_abline(data = ablines, aes(slope = a, intercept = b, color = Group), alpha = 0.01)

            print(p_Regression)
        })
    })
## run app
shinyApp(ui = ui, server = server)
