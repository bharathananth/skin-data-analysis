#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(magrittr)
library(tidyr)
library(dplyr)


#load data
results <- readRDS("data/results.rds")
y <- readRDS("data/rawdata.rds")
experiment <- readRDS("data/experiment.rds")
print(dim(y))

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Human skin dataset"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        
         sliderInput("cutoff", "FDR cutoff", 0.001, 0.25,0.05, step=0.025),
         
         uiOutput("transcriptControls")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        plotOutput("timecourse")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$transcriptControls <- renderUI({
     circ_transcripts <- filter(results, adj_P_Val < input$cutoff)
     tagList(
     selectInput("transcript",
                 "Transcript:",
                 circ_transcripts["ProbeName"])
     )
   })
   
   output$timecourse <- renderPlot({
     # chose transcript base input$transcript from ui.R
     data <- data.frame(value = as.numeric(y$E[input$transcript,]),
                        time = experiment$time,
                        tissue = experiment$tissue,
                        subject = experiment$subject)
     
     title <- filter(results, ProbeName == input$transcript) %$% Symbol
     
     ggplot(data, aes(x=time, y = value, color=subject)) + geom_line() + geom_point() + 
       facet_wrap(~tissue, scales = "free_y", ncol = 2) +
       theme_bw() + theme(aspect.ratio = 1) + ggtitle(title)
     
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

