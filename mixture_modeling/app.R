library(rhandsontable)
library(rjags)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(tidybayes)
library(shiny)
theme_set(theme_cowplot(font_size = 18))

ui <- shinyUI(fluidPage(
    title = "Gaussian mixture",
    titlePanel("Gaussian mixture"),
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
                                 plotOutput("plotOver",height = "300px"),
                                 helpText("The last two plots can help evaluate if the variability in the data can be attributed to other nuisance factors."),
                                 plotOutput("plotPost",height = "400px"),
                                 helpText("The second plot shows the credible probability distributions of the difference between the group means without the nuisance factors. This information can be used to determine how significant the difference is, but also, to evaluate the magnitude of the difference."),
                                 plotOutput("plotDiffGroup",height = "300px")
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
                Data = round(c(rnorm(20,10,1.5),rnorm(5,12,0.5)),digits = 1),
                Group = c(rep("A",20), rep("B",5)),
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
            DF$Group = as.factor(DF$Group)
            p <- ggplot(DF) + geom_dotplot(aes(x = Data, fill = Group), stackgroups = TRUE, binwidth = 0.5, method = "histodot") + expand_limits(x = c(min(DF$Data)-1, max(DF$Data)+1)) + scale_y_continuous(limits = c(0,0.4))
            print(p)
        })

        samples <- reactive({
            withProgress(message = 'Calculation in progress...', value = 0, {
                DF <- values[["DF"]]

                Data = DF$Data
                Group = as.factor(DF$Group)
                N = length(DF$Data)
                clust = as.numeric(Group)
                # clust[which.min(Data)] = 1 # smallest value assigned to cluster 1
                # clust[which.max(Data)] = 2 # highest value assigned to cluster 2

                data = list(Data = DF$Data,
                            N = N,
                            Nclust = 2,
                            meanData = mean(Data),
                            clust = clust,
                            onesRepNclust = rep(1,2))

                model <-
                    jags.model(
                        "mixture_modeling.bug",
                        data = data,
                        n.chains = 4
                    )

                update(model, n.iter = 10000)

                samples <-
                    coda.samples(
                        model,
                        variable.names = c("muOfClust", "tauOfClust", "pClust"),
                        n.iter = 5000,
                        thin = 100
                    )

                samples = recover_types(samples, DF)
            })
        })

        output$plotOver = renderPlot({
            samples = samples()
            DF <- values[["DF"]]
            DF$Group = as.factor(DF$Group)
            lc = c("#F8766D","#00BFC4")
            draws = samples %>%
                spread_draws(muOfClust[clust],tauOfClust[clust],pClust[clust])
            p <- ggplot(DF) + geom_dotplot(aes(x = Data, fill = Group), stackgroups = TRUE, binwidth = 0.5, method = "histodot") + expand_limits(x = c(min(DF$Data)-1, max(DF$Data)+1)) + scale_y_continuous(limits = c(0,0.4))
            x = seq(min(DF$Data) - 1, max(DF$Data) + 1, 0.1)
            for (i in seq(1,length(draws$clust),by=5)) {
                d_data = data.frame(x = x, y = draws$pClust[i] * dnorm(x,draws$muOfClust[i], sqrt(1/draws$tauOfClust[i])))
                p <- p + geom_line(data = d_data, aes(x = x, y = y), color = lc[draws$clust[i]], alpha = 0.1, size = 1)
            }
            print(p)
        })

        output$plotDiffGroup = renderPlot({
            samples = samples()
            p_diff_Group_Mean = samples %>%
                spread_draws(muOfClust[Group]) %>%
                compare_levels(muOfClust, by = Group) %>%
                ggplot(aes(y = fct_rev(Group), x = muOfClust)) + geom_eyeh(fill = 'firebrick') + geom_vline(xintercept = 0, linetype = 'dashed') + ylab('Group')
            p_diff_Group_Tau = samples %>%
                spread_draws(tauOfClust[Group]) %>%
                compare_levels(tauOfClust, by = Group) %>%
                ggplot(aes(y = fct_rev(Group), x = sqrt(1/tauOfClust))) + geom_eyeh(fill = 'firebrick') + geom_vline(xintercept = 0, linetype = 'dashed') + ylab('Group')
            p_diff_Group_Prop = samples %>%
                spread_draws(pClust[Group]) %>%
                compare_levels(pClust, by = Group) %>%
                ggplot(aes(y = fct_rev(Group), x = pClust)) + geom_eyeh(fill = 'firebrick') + geom_vline(xintercept = 0, linetype = 'dashed') + ylab('Group')
            p = plot_grid(p_diff_Group_Mean, p_diff_Group_Tau, p_diff_Group_Prop, ncol = 3)
            print(p)
        })
        output$plotPost = renderPlot({
            samples = samples()
            p_Mean = samples %>%
                spread_draws(muOfClust[Group]) %>%
                ggplot(aes(y = fct_rev(Group), x = muOfClust)) + geom_eyeh(fill = 'firebrick') + xlab('Mean')

            p_Std = samples %>%
                spread_draws(tauOfClust[Group]) %>%
                ggplot(aes(y = fct_rev(Group), x = sqrt(1/tauOfClust))) + geom_eyeh(fill = 'firebrick') + xlab('Std')

            p_Prop = samples %>%
                spread_draws(pClust[Group]) %>%
                ggplot(aes(y = fct_rev(Group), x = pClust)) + geom_eyeh(fill = 'firebrick') + xlab('Prop')

            p = plot_grid(p_Mean, p_Std, p_Prop, ncol = 3)
            print(p)
        })
    })
## run app
shinyApp(ui = ui, server = server)
