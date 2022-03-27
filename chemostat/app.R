library(shiny)
library(deSolve)
library(ggplot2)
library(cowplot)
library(reshape2)
theme_set(theme_cowplot())

ui <- fluidPage(
  titlePanel("Chemostat"),
  sidebarLayout(
    sidebarPanel(
      h4("Biological parameters"),
      sliderInput('q_s_max', 'Max. substrate consumption rate (Qsmax)', min = 0, max = 2, value = 1, step = 0.1, round = FALSE),
      sliderInput('m_s', 'Metabolic maintenance cost (Ms)', min = 0, max = 0.2, value = 0.1, step = 0.01, round = FALSE),
      sliderInput('k_s', 'Inverse substrate affinity (1/ks)', min = 0, max = 1, value = 0.5, step = 0.05, round = FALSE),
      sliderInput('Y_xs_max', 'Maximum biomass yield (Ymax)', min = 0, max = 1, value = 0.5, step = 0.05, round = FALSE),
      hr(),
      h4("Reactor parameters"),
      sliderInput('C_s_0', 'Starting substrate concentration ([S] t=0)', min = 0, max = 20, value = 10, step = 1, round = FALSE),
      sliderInput('C_x_0', 'Starting cell concentration ([C] t=0)', min = 0, max = 20, value = 10, step = 1, round = FALSE),
      sliderInput('D', 'Dilution rate (D)', min = 0, max = 0.25, value = 0.125, step = 0.0125, round = FALSE),
      sliderInput('C_s_in', 'Limiting substrtate concentration ([S]in)', min = 0, max = 20, value = 10, step = 1, round = FALSE),
      sliderInput('cost', 'Operating cost', min = 0, max = 10, value = 5, step = 0.05, round = FALSE)
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Model",
                           h4("Diagram"),
                           img(src = "chemostat.png", height = '500px'),
                           h4("Model equations"),
                           img(src = "equations.png", height = '300px'),
                           h4("Model parameters"),
                           p("[S]in: input substrate concentration"),
                           p("[S]: substrate concentraion in reactor"),
                           p("[C]: cell concentration in reactor"),
                           p("D: dilution rate"),
                           p("mu: cell growth rate"),
                           p("Ymax: maximum cell yield"),
                           p("Ms: metabolic maintenence cost"),
                           p("Qsmax: maximum consumption rate"),
                           p("Qs: consumption rate")
                  ),
                  tabPanel("Dynamics",
                           h4("Chemostat dynamics"),
                           plotOutput("plot_output", height = "300px"),
                           plotOutput('plot_yield', height = "300px")),
                  tabPanel("Phase plane",
                           h4("Chemostat phase plane"),
                           plotOutput('plot_phase', height = "600px")),
                  tabPanel("Optimal performance", h4("Optimal chemostat performance"),
                           plotOutput('plot_heatmaps', height = "800px"))
      )
    )
  )
)



server <- function(input, output, session) {

  session$onSessionEnded(stopApp)
  times <- seq(0, 120, by = 0.01)

  out <- reactive({
    parameters <- c(D = input$D, C_s_in = input$C_s_in, q_s_max = input$q_s_max, m_s = input$m_s, k_s = input$k_s, Y_xs_max = input$Y_xs_max)
    state <- c(C_s = input$C_s_0, C_x = input$C_x_0)
    chemostat <- function(t, state, parameters) {
      with(as.list(c(state, parameters)),{
        q_s <- q_s_max * (C_s / (k_s + C_s))
        mu <- (q_s - m_s) * Y_xs_max
        dC_s <- D * C_s_in - D * C_s - q_s * C_x
        dC_x <- C_x * mu - D * C_x
        list(c(dC_s, dC_x))
      })
    }

    out <- as.data.frame(ode(y = state, times = times, func = chemostat, parms = parameters))
  })

  output$plot_yield <- renderPlot({
    data = out()
    q_s <- input$q_s_max * (data$C_s / (input$k_s + data$C_s))
    Y <- (q_s - input$m_s) * input$Y_xs_max / q_s
    data = melt(cbind(data,q_s,Y), id.vars = "time", measure.vars = c("q_s","Y"))
    p <- ggplot(data) + geom_line(aes(x = time, y = value, color = variable), size = 2) + expand_limits(y = 0) + xlim(0, 60) + scale_color_brewer(name = "Variables", labels = c("Substrate consumption rate", "Yield"), palette = "Set2") + theme(legend.position = "top")
    print(p)
  })

  output$plot_output <- renderPlot({
    data = melt(data = out(), id.vars = "time", measure.vars = c("C_s", "C_x"))
    p <- ggplot(data) + geom_line(aes(x = time, y = value, color = variable), size = 2) +
      expand_limits(y = 0) + xlim(0, 60) +
      scale_color_brewer(name = "Variables", labels = c("Substrate concentration", "Cells concentration"), palette = "Set1") + theme(legend.position = "top")
    print(p)
  })

  output$plot_phase <- renderPlot({

    data = out()

    C_s_seq = seq(0,20,0.02)
    C_x_seq = seq(0,20,0.02)

    C_x_eq = -input$D * (C_s_seq - input$C_s_in) / (input$m_s + input$D / input$Y_xs_max)
    C_s_eq = (sqrt(input$D^2 * (input$k_s + input$C_s_in)^2 + 2 * input$D * (input$k_s - input$C_s_in) * input$q_s_max * C_x_seq + input$q_s_max^2 * C_x_seq^2) - input$D * (input$k_s - input$C_s_in) - input$q_s_max * C_x_seq) / (2 * input$D)

    ss = as.data.frame(cbind(C_s_seq,C_x_seq,C_s_eq,C_x_eq))

    grid = as.data.frame(expand.grid(seq(0, 20), seq(0, 30)))
    colnames(grid)  <- c("C_s", "C_x")

    q_s <- input$q_s_max * (grid$C_s / (input$k_s + grid$C_s))
    mu <- (q_s - input$m_s) * input$Y_xs_max
    dC_s <- (input$D * input$C_s_in - input$D * grid$C_s - q_s * grid$C_x)
    dC_x <- (grid$C_x * mu - input$D * grid$C_x)

    dM = sqrt(dC_s^2 + dC_x^2)

    dC_s = dC_s/dM/2
    dC_x = dC_x/dM/2

    dM = dM/max(dM)

    grid = cbind(grid, dC_s, dC_x, dM)

    p <- ggplot() +
      geom_segment(data = grid, aes(x = C_s, y = C_x, xend = C_s+dC_s, yend = C_x+dC_x), alpha = dM, size = 1) +
      geom_point(data = grid, aes(x = C_s+dC_s, y = C_x+dC_x), alpha = dM, size = 1) +
      geom_line(data = ss, aes(x = C_s_seq, y = C_x_eq), size = 0.5) +
      geom_line(data = ss, aes(x = C_s_eq, y = C_x_seq), size = 0.5) +
      geom_point(data = out(), aes(x = C_s, y = C_x, color = log(time)), size = 3) +
      xlim(0, 20) + ylim(0, max(20, out()$C_x)) + theme(legend.position = "none", aspect.ratio = 1) + xlab("Subtrate conc.") + ylab("Cells conc.")
    print(p)
  })

  output$plot_heatmaps <- renderPlot({

    D = seq(0,0.25,0.01)
    C_s_in = seq(0,20,0.8)
    data = expand.grid(D = D,C_s_in = C_s_in)

    q_s = data$D / input$Y_xs_max + input$m_s
    C_s = input$k_s * q_s / (input$q_s_max - q_s)
    C_x = data$D / q_s * (data$C_s_in - C_s)

    C_x[q_s > input$q_s_max * (C_s/(C_s + input$K_s))] = 0
    C_x[C_s < 0] = 0
    C_x[q_s <= 0] = 0
    C_x[C_s > data$C_s_in] = 0
    C_x[is.nan(C_s)] = 0

    Y = C_x / data$C_s_in
    Y[!is.finite(Y)] = 0

    P = C_x * D / (D * data$C_s_in^1.5 + input$cost)
    P[D <= 0] = 0

    data  = cbind(data,Y,P)

    pY <- ggplot(data, aes(x = D, y = C_s_in)) + geom_tile(aes(fill = Y)) + scale_fill_distiller(name = "Yield", palette = "Blues") + theme(aspect.ratio = 1) + xlab("Dilution rate") + ylab("Limiting substrtate conc.")
    pP <- ggplot(data, aes(x = D, y = C_s_in)) + geom_tile(aes(fill = P)) + scale_fill_distiller(name = "Return", palette = "Greens")  + theme(aspect.ratio = 1) + xlab("Dilution rate") + ylab("Limiting substrtate conc.")

    p  = plot_grid(pY,pP,ncol =1)

    #p <- ggplot(data) + geom_line(aes(x = time, y = value, color = variable), size = 2) + expand_limits(x = 0, y = 0) + scale_color_brewer(name = "Variables", labels = c("Substrate concentration", "Cells concentration"), palette = "Set1") + theme(legend.position = "top")
    print(p)
  })
}

shinyApp(ui = ui, server = server)


