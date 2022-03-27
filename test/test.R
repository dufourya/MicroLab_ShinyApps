library(rhandsontable)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(shiny)

ui <- shinyUI(
  fluidPage(
    titlePanel("Models"),
    sidebarLayout(
      sidebarPanel(
        helpText("Choose model"),
        wellPanel(
          h3("Options"),
          radioButtons("model", "Model", c("Average", "Independent Fits", "Hierarchical fit"))
        ),
        br(), 
        wellPanel(
          h3("Fit"), 
          actionButton("fit", "Fit")
        )        
      ),
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Data table", rHandsontableOutput("hot")),
                    tabPanel("Plot", plotOutput("plot")),
                    tabPanel("Summary", verbatimTextOutput("summary"))
        )
      )
    )
  )
)

cell_growth <- function(time_points, initial_OD, proportion_dead_cells, carrying_capacity,growth_rate){
  
  OD = initial_OD * proportion_dead_cells + carrying_capacity * 2^(time_points / growth_rate) / (
    2^(time_points/growth_rate) - 1 + carrying_capacity / (initial_OD * (1 - proportion_dead_cells)))
  
  return(OD)
}

time_points = seq(1,500,by = 30)
initial_OD =  0.01
proportion_dead_cells = 0.25
carrying_capacity = 1.2
growth_rate = 30

DF <- data.frame(Time = NULL, OD = NULL, Replicate = NULL, stringsAsFactors = FALSE)

for (i in 1:3) {
  t = time_points + (runif(length(time_points)) * 20)
  OD = cell_growth(t, initial_OD * 2*(runif(1) + 0.5), proportion_dead_cells * 2*(runif(1) + 0.5), carrying_capacity + 0.1*runif(1), growth_rate + 2*runif(1))
  DF = rbind(DF, data.frame(Time = t, OD = OD, Replicate = rep(i,length(OD))))
}

server <- shinyServer(function(input, output) {
  
  values <- reactiveValues()
  
  ## Handsontable
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
      rhandsontable(DF, useTypes = as.logical(TRUE), stretchH = "all", readOnly = FALSE)
  })
  
  output$plot <- renderPlot({
    DF <- values[["DF"]]
    DF$Replicate = as.factor(DF$Replicate)
    p <- ggplot(DF) + geom_point(aes(x = Time, y = OD, color = Replicate), size = 4)
    print(p)
  }
  )
  
  ## Fit
  observeEvent(input$fit, {
    finalDF <- isolate(values[["DF"]])
  })
  
})

## run app 
runApp(list(ui = ui, server = server))

