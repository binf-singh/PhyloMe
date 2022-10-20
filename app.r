# Load libraries
library("shiny")
library("shinyWidgets")  
library("ape")
library("seqinr")
library("Biostrings") 
library("tidyverse") 
library("phangorn") 
library("phytools")
library("msa") 
library("stringr")
library("taxize")
library("ggtree") 
library("ggplot2")
        

# Retrieve UI and Server files
source("ui.R")
source("server.R")

# Create shiny app
shinyApp(ui = ui, server = server)
