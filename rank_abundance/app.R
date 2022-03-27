library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

ui <- fluidPage(
  titlePanel("Quantifying Microbial Community Diversity"),
  sidebarLayout(
    sidebarPanel(
      h3("Parameters"),
      sliderInput('num_OTU', 'Number of OTUs', min = 1, max = 50, value = 25, step = 1, round = TRUE),
      sliderInput('evenness', 'OTU abundance evenness', min = 0, max = 1, value = 0.5, step = 0.1, round = FALSE),
      sliderInput('seq_depth', 'Sequencing depth (log10)', min = 1, max = 5, value = 3, step = 0.1, round = FALSE),
      hr(),
      actionButton("save_button", "Save plot")
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Introduction",
                           p("You learned in class that microbes are everywhere and scientists are interested in uncovering their functional capabilities and the potential roles they are playing in their environments. You also learned that 99% of microbes cannot be cultivated in the lab, so we must look to alternative methods, such as high-throughput sequencing to study these microbes."),
                           p("But how do we know that when conducting these studies/surveys, we are capturing as much of the diversity of microbes present as possible? How do we know that we recovered the majority of microbial members present in our system?"),
                           p("Scientists use two methods to address those questions: rarefaction curves and rank-abundance curves."),
                           p("Rarefaction curves enumerate the number of unique operational taxinomic units (OTUs) found in the sequences. Sequencing depth  is the number of usable reads per sample obtained from sequencing. If you sequenced deeply, you have many genomic reads per sample."),
                           p("Rank-abundance curves show the relative abundances (i.e. proportions) of OTUs in your sample, ordered from most to least abundant. These curves can tell you both about the richness (# unique OTUs) and evenness (how 'evenly' distributed those OTUs are) of your community."),
                           p("In this interactive exercise, we are going to explore how community complexity (# OTUs in your community) and community evenness interact with sequencing effort to determine whether we have comprehensively sampled our community."),
                           p("As you change the values of Number of OTUs, OTU Abundance Evenness, and Sequencing Depth, think about these questions:"),
                           p("- How do rarefaction curves and rank-abundance plots change as communities become more complex?"),
                           p("- Why does having a community with uneven OTU abundances make it harder to capture the entire diversity of the system?"),
                           p("- Is aiming for the highest sequencing depth always the best way to go?"),
                           p("- What do rarefaction curves with long right-tails say about your community?"),
                           p("Note: the Save plot button will save that plot and layer it on future plots so that you can compare the effects of parameters.")),
                  tabPanel("Plots",
                           h3("Rarefaction"),
                           plotOutput("plot_rarefaction", height = "300px"),
                           h3("Rank-abundance"),
                           plotOutput("plot_rank", height = "300px"))
      )
    )
  )
)

p = ggplot()

server <- function(input, output, session) {
  session$onSessionEnded(stopApp)

  OTU = reactive({
    OTU_freq = exp(seq(-10,0,0.1))^(1-input$evenness)
    num_OTU = input$num_OTU

    OTU_freq = sample(OTU_freq, num_OTU)
    OTU_freq = OTU_freq/sum(OTU_freq)

    OTU = sample(seq(1,num_OTU), 10^5, replace = TRUE, prob = OTU_freq)
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
      expand_limits(y = 0, x = 10^5) + scale_x_continuous(trans = 'log10', name = 'Sequencing depth (log10(# genomic reads))') +
      scale_y_continuous(name = 'Number of unique OTUs')
    print(p)

  })

  saved_data <- eventReactive(input$save_button, {
    list(rare = rare(), table = OTU_table())
  })

}

shinyApp(ui = ui, server = server)
