library(shiny)
library(bslib)
library(tinyplot)
source("calculator.R")

ui <- page_navbar(
  title = "Powered for Paediatric Safety",
  nav_panel(
    "Information",
    "Some text describing the model, inputs and calculations"
  ),
  nav_panel(
    "Calculator",
    layout_columns(
      card(
        layout_columns(
          card(
            card_header("Reference Population"),
            numericInput("ref_cx", "Control rate", value = 0.01, min = 0, max = 1, step = 0.01),
            numericInput("ref_tx", "Treated rate", value = 0.03, min = 0, max = 1, step = 0.01)
          ),
          card(
            card_header("Prior Distributions"),
              ("Beta Prior for Control"),
              layout_columns(
                shinyWidgets::numericInputIcon("prior_cx_a", label = NULL, icon= "\u03b1", value = 1, min = 0, step = 0.5),
                shinyWidgets::numericInputIcon("prior_cx_b", label = NULL, icon= "\u03b2", value = 1, min = 0, step = 0.5)
              ),
            ("Beta Prior for Treated"),
            layout_columns(
              shinyWidgets::numericInputIcon("prior_tx_a", label = NULL, icon= "\u03b1", value = 1, min = 0, step = 0.5),
              shinyWidgets::numericInputIcon("prior_tx_b", label = NULL, icon= "\u03b1", value = 1, min = 0, step = 0.5)
            )
          ),
          card(
            card_header("Safety Database"),
            shinyWidgets::numericRangeInput("db_cx", "Control Database (#Events / N)", value = c(0,0), separator = "/"),
            shinyWidgets::numericRangeInput("db_tx", "Treated Database (#Events / N)", value = c(0,0), separator = "/"),
            # numericInput("ped_cx_n", "Control N", value = 0, min = 0),
            # numericInput("ped_cx_event", "Control Events", value = 0, min = 0),
            # numericInput("ped_tx_n", "Treated N", value = 0, min = 0),
            # numericInput("ped_tx_event", "Treated Events", value = 0, min = 0)
          ),
          card(
            card_header("New Study"),
            numericInput("new_n", "Sample Size", value = 40, min = 0),
            shinyWidgets::numericRangeInput("rand", "Randomization", separator = ":", value = c(1, 1)),
            numericInput("fold_increase", "Fold increase over reference difference", value = 2, min = 1)
          )
        )
      ),
      card(
        layout_columns(
          plotOutput("curve"),
          tableOutput("table"),
          col_widths = c(6, 6)
        )
      ),
      col_widths = c(12, 12)
    )
  )
)

server <- function(input, output) {

  # Constrain inputs dynamically
  # observeEvent(input$ped_cx_n, {
  #   updateSliderInput(inputId = "ped_cx_event", max = input$ped_cx_n)
  # })
  # observeEvent(input$ped_tx_n, {
  #   updateSliderInput(inputId = "ped_tx_event", max = input$ped_tx_n)
  # })

  prob_table <- reactive({
    prob_rule_out2(
      ref_tx = input$ref_tx,
      ref_cx = input$ref_cx,
      fold_increase = input$fold_increase,
      prior_tx_a = input$prior_tx_a,
      prior_tx_b = input$prior_tx_b,
      prior_cx_a = input$prior_cx_a,
      prior_cx_b = input$prior_cx_b,
      ped_cx_event = input$db_cx[1],
      ped_cx_n = input$db_cx[2],
      ped_tx_event = input$db_tx[1],
      ped_tx_n = input$db_tx[2],
      new_n = input$new_n,
      rand = input$rand
    )
  })

  output$table <- renderTable({
    prob_table()
  }, digits = 4)

  output$curve <- renderPlot({
    df <- prob_table()

    # plot(
    #   x = df$N_t + df$N_c,
    #   y = df$Probability,
    #   type = "n",
    #   ylim = c(0, 1),
    #   xlab = "Study Size",
    #   ylab = "Probability"
    # )
    #
    # idf <- split(df, df$`Incidence Difference Factor`)
    # for (i in seq_along(idf)) {
    #
    #   print(points(
    #     x = idf[i]$N_t + idf[i]$N_c,
    #     y = idf[i]$Probability,
    #     type = "b"
    #   ))
    # }

    tinyplot(
      x = df$N_t + df$N_c,
      y = df$Probability,
      by = df$`Incidence Difference Factor`,
      type = "b",
      # ylim = c(0, 1),
      xlab = "Study Size",
      ylab = "Probability"
    )
  })

}

shinyApp(ui = ui, server = server)
