library(shiny)
library(bslib)
source("calculator.R")

ui <- page_fillable(
  title = "Powered for Paediatric Safety",
  layout_columns(
    card(
      card_header("Inputs"),
      layout_columns(
        card(
          card_header("Reference Population"),
          numericInput("ref_cx", "Control rate", value = 0.01, min = 0, max = 1, step = 0.01),
          numericInput("ref_tx", "Treated rate", value = 0.03, min = 0, max = 1, step = 0.01)
          ),
        card(
          card_header("Prior Information"),
          numericInput("prior_cx_a", "Control alpha", value = 1, min = 0, step = 0.5),
          numericInput("prior_cx_b", "Control beta", value = 1, min = 0, step = 0.5),
          numericInput("prior_tx_a", "Treated alpha", value = 1, min = 0, step = 0.5),
          numericInput("prior_tx_b", "Treated beta", value = 1, min = 0, step = 0.5)
        ),
        card(
          card_header("Safety Database"),
          numericInput("ped_cx_n", "Control N", value = 0, min = 0),
          numericInput("ped_cx_event", "Control Events", value = 0, min = 0),
          numericInput("ped_tx_n", "Treated N", value = 0, min = 0),
          numericInput("ped_tx_event", "Treated Events", value = 0, min = 0)
        ),
        card(
          card_header("New Study"),
          numericInput("ped_n_cx", "N control arm", value = 20, min = 0),
          numericInput("ped_n_tx", "N treated arm", value = 20, min = 0),
          numericInput("fold_increase", "Fold increase over reference difference", value = 2, min = 1)
        )
      )
    ),
    card(
      card_header("Results"),
      tableOutput("table")
    ),
    col_widths = c(12, 12)
  )
)

server <- function(input, output) {

  # Constrain inputs dynamically
  observeEvent(input$ped_cx_n, {
    updateSliderInput(inputId = "ped_cx_event", max = input$ped_cx_n)
  })
  observeEvent(input$ped_tx_n, {
    updateSliderInput(inputId = "ped_tx_event", max = input$ped_tx_n)
  })


  output$table <- renderTable({
    # browser()
    prob_rule_out2(
      ref_tx = input$ref_tx,
      ref_cx = input$ref_cx,
      fold_increase = input$fold_increase,
      ped_n_cx = input$ped_n_cx,
      ped_n_tx = input$ped_n_tx,
      prior_tx_a = input$prior_tx_a,
      prior_tx_b = input$prior_tx_b,
      prior_cx_a = input$prior_cx_a,
      prior_cx_b = input$prior_cx_b,
      ped_cx_event = input$ped_cx_event,
      ped_cx_n = input$ped_cx_n,
      ped_tx_event = input$ped_tx_event,
      ped_tx_n = input$ped_tx_n
    )
  }
  )
}

shinyApp(ui = ui, server = server)


#
# ref_tx = 0.06,
# ref_cx = 0.03,
# fold_increase = 3,
# rand_ratio = 1/1,
# ped_n = 200,
#  = 1,
#  = 1,
#  = 1,
# prior_cx_b = 1,
#  = 0,
#  = 0,
#  = 0,
#  = 0
