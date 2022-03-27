library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

ui <- fluidPage(
  titlePanel("Differential gene expression in RNA-Seq Analysis"),
  sidebarLayout(
    sidebarPanel(
      h3("Parameters"),
      sliderInput('gene_evenness', 'Gene expression evenness', min = 0, max = 1, value = 0.5, step = 0.1, round = FALSE),
      sliderInput('fold_change', 'Fold change (condition A vs. B)', min = 0, max = 1, value = 0.5, step = 0.1, round = FALSE),
      sliderInput('seq_depth', 'Sequencing depth (log10)', min = 1, max = 6, value = 3, step = 0.1, round = FALSE),
      sliderInput('type_I', 'Type I error (p-value)', min = 0, max = 0.09, value = 0.02, step = 0.005, round = FALSE),
      sliderInput('type_II', 'FDR stringency', min = 0, max = 0.06, value = 0.03, step = 0.005, round = FALSE),
      hr(),
      actionButton("save_button", "Save plot")
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Introduction",
                  p("RNA-Seq can be a powerful tool to quantify gene expression levels and determine if they vary across tissue types, conditions, and between groups. We conduct statistics to determine precisely which genes are under or over-expressed in certain groups."),
                  p("Because RNA-Seq studies commonly include thousands of genes, this means we are conducting thousands of tests simultaneously, and the likelihood that we observe a statistically significant result simply due to chance increases. We must control for this Type 1 error (incorrect rejection of a true null hypothesis) which is also known as α or p-value. In genomics, we can apply a false discovery rate (FDR) correction to our p-values, and can set this correction to a desired level of stringency. Being too lenient or too stringent is going to affect both the rate of false positives and our power to detect differences."),
                  p("Another factor that will impact our findings and detection of differentially expressed genes is the evenness in the expression of genes themselves: do you have many highly-expressed genes and many low-expressed genes? Or genes that are found at more equal abundances? You also have to think about the magnitude of difference (fold change) you expect among your treatments/groups. Do you expect the differences in gene expression among conditions to be high or low? Lastly, as you learned in the prior Shiny App, sequencing coverage/depth will be important; you want your coverage to be high enough so that you have sequences that cover the entire length of your transcript and no areas are missed."),
                  p("In this App you will get to explore how the evenness of gene expression, fold change between conditions, sequencing depth, type 1 error, and FDR levels affect findings regarding differential gene expressions between two hypothetical treatment conditions (A vs. B).")),
                  tabPanel("Plots",
                           h3("MA-Plot"),
                           plotOutput("plot_ma", height = "300px"),
                           h3("Volcano Plot"),
                           plotOutput("plot_volcano", height = "300px"))
      )
    )
  )
)

##YANN'S CODE FROM LAST APP TO HELP ME BUILD CODE FOR NEW APP
p = ggplot()

server <- function(input, output, session) {
  session$onSessionEnded(stopApp)

  OTU = reactive({
    OTU_freq = exp(seq(-10,0,0.1))^(1-input$evenness)
    num_OTU = input$num_OTU

    OTU_freq = sample(OTU_freq, num_OTU)
    OTU_freq = OTU_freq/sum(OTU_freq)

    OTU = sample(seq(1,num_OTU), 10^6, replace = TRUE, prob = OTU_freq)
  })

  OTU_sub = reactive({
    sample(OTU(),10^input$seq_depth,replace = FALSE)
  })

  rare = reactive({
    OTU = OTU_sub()
    rare = tibble(seq_num = 10^(seq(1,input$seq_depth,0.1)), otu_num = 0)
    for (i in seq(1,length(rare$seq_num))) {
      for (j in 1:100){
        rare$otu_num[i] = rare$otu_num[i] + length(unique(sample(OTU,i,replace = FALSE)))
      }
    }
    rare$otu_num = rare$otu_num/100
    rare

  })

  OTU_table = reactive({
    OTU = OTU_sub()
    OTU_table = as.data.frame(table(OTU))
    OTU_table = OTU_table[order(OTU_table$Freq,decreasing = TRUE),]
    OTU_table$Freq = OTU_table$Freq/sum(OTU_table$Freq)*100
    OTU_table$OTU = factor(OTU_table$OTU, levels = OTU_table$OTU)
    OTU_table
  })

  output$plot_rank <- renderPlot({

    OTU_table = OTU_table()
    num_OTU = input$num_OTU
    if (input$save_button > 0) {
    saved_data = saved_data()
    saved_table = saved_data$table
    p = p + geom_step(data = saved_table, aes(x = seq(1,length(saved_table$Freq)), y = Freq),direction='vh',color='firebrick')
    }
    p = p + geom_step(data = OTU_table, aes(x = seq(1,length(OTU_table$Freq)), y = Freq),direction='vh') +
      expand_limits(y = 0, x = num_OTU) +
      scale_x_continuous(breaks = seq(1,length(OTU_table$Freq)), labels = OTU_table$OTU, name = 'OTU ID') +
      scale_y_continuous(name = 'Relative percentage')
    print(p)

  })

  output$plot_rarefaction <- renderPlot({

    num_OTU = input$num_OTU
    rare = rare()
    if (input$save_button > 0) {
    saved_data = saved_data()
    saved_rare = saved_data$rare
    p = p + geom_point(data = saved_rare, aes(x = seq_num, y = otu_num), color = 'firebrick')
    }
    p = p + geom_point(data = rare, aes(x = seq_num, y = otu_num)) +
      geom_hline(yintercept = num_OTU, linetype = 'dotted') +
      expand_limits(y = 0, x = 10^6) + scale_x_continuous(trans = 'log10', name = 'Sequencing depth') +
      scale_y_continuous(name = 'Number of unique OTUs')
    print(p)

  })

  saved_data <- eventReactive(input$save_button, {
    list(rare = rare(), table = OTU_table())
  })

}

shinyApp(ui = ui, server = server)
