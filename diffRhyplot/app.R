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
suppressPackageStartupMessages(library(SummarizedExperiment))
library(gridExtra)
library(limma)
library(edgeR)


#load data
se <- readRDS("../data/GSE52333.RDS")
experiment_1 <- colData(se)[, c("harvest.timepoint.ch1", "diet.ch1")]
experiment_1 <- as.data.frame(experiment_1)
experiment_1$time <- as.numeric(sub("ZT", "", experiment_1$harvest.timepoint.ch1))
experiment_1$group <- factor(make.names(experiment_1$diet.ch1))
y_1 <- assay(se, "expr")

se <- readRDS("../data/PRJNA428303.RDS")
experiment_2 <- colData(se)[, c("time", "group")]
experiment_2 <- as.data.frame(experiment_2)
experiment_2$time <- as.numeric(sub("ZT", "", experiment_2$time))
experiment_2$group <- factor(experiment_2$group)
y_2 <- assay(se, "countsFromAbundance")


# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("compareRhythms datasets"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         textInput("gene", "ENSEMBL ID", value="ENSMUSG00000020893")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        plotOutput("timecourse")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$timecourse <- renderPlot({
     # chose transcript base input$transcript from ui.R
     data <- data.frame(value = y_1[input$gene,],
                        time = experiment_1$time,
                        condition = experiment_1$group)
     
     mean_data <- data %>% group_by(condition, time) %>% summarize(value=mean(value))
     
     f1 <- ggplot(data, aes(x=time, y = value, color=condition)) + geom_point() + 
       geom_line(data=mean_data) + 
       theme_bw() + theme(aspect.ratio = 1)
     
     data <- data.frame(value = y_2[input$gene,],
                        time = experiment_2$time,
                        condition = experiment_2$group)
     
     mean_data <- data %>% group_by(condition, time) %>% summarize(value=mean(value))
     
     f2 <- ggplot(data, aes(x=time, y = value, color=condition)) + geom_point() + 
        geom_line(data=mean_data) + 
        theme_bw() + theme(aspect.ratio = 1)
     
     grid.arrange(f1, f2, ncol=1)
     
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

