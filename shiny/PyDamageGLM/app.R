#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

mod = readRDS("data/pydamage_glm_model.rds")

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("PyDamage Accuracy prediction"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            numericInput("cov",
                        "Contig Mean coverage",
                        min = 0.0,
                        max = 99999999999.0,
                        value = 25.3,
                        step=0.01),
            numericInput("length",
                         "Contig Lenth",
                         min = 0.0,
                         max = 99999999999,
                         value = 10000,
                         step=1),
            numericInput("damage",
                         "Damage on 5' end",
                         min = 0.0,
                         max = 0.99,
                         value = 0.3,
                         step=0.01)
        ),

        # Show a plot of the generated distribution
        mainPanel(
            htmlOutput("prediction")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$prediction <- renderText({
        input_data = data.frame(actual_cov=input$cov, damage=input$damage, contiglength=input$length)
        prediction = round(predict(mod, input_data, type='response'), 3)
        paste("<b>Predicted accuracy:</b>", prediction)
    })

}
# Run the application 
shinyApp(ui = ui, server = server)
