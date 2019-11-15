# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.

# Author: Julien Bryois
# Date: 14.11.2019

library(shiny)
library(tidyverse)
library(broom)

# Load file with helper functions
source(file = "functions.R")

#options to increase max size of the file to be loaded
options(shiny.maxRequestSize=90*1024^2) 

# Define server logic
shinyServer(function(input, output) {

    # Load list of genes of interest (user input)
    gene_list <- reactive({
        inFile <- input$file1
        req(inFile)
        if(!input$sep%in%c('xlsx','xls')){
            tbl <- read.csv(inFile$datapath, header=input$header, sep=input$sep,comment = "#",stringsAsFactors = FALSE)
        } 
        if(input$sep=='xlsx'){
            tbl <- as.data.frame(read_xlsx(inFile$datapath,col_names=input$header))
        }
        if(input$sep=='xls'){
            tbl <- as.data.frame(read_xls(inFile$datapath,col_names=input$header))
        }
        return(tbl)
    }) 
    
    # Load selected dataset (e.g. GTEx v8)
    dataset <- reactive({
        inFile <- input$select
        req(inFile)
        d <- read_tsv(inFile)
    }) 
    
    # Get regression results
    compute_reg <- reactive({
        reg(dataset(),gene_list())
    }) 
    
    # Get fisher results
    compute_fisher <- reactive({
        fisher_test(dataset(),gene_list())
    }) 
    
    # Get EWCE results
    compute_EWCE <- reactive({
        ewce(dataset(),gene_list())
    }) 
    
    # Plot regression results
    output$reg <- renderPlot({
        plot_results(compute_reg(),input$plot_type)
    }, height = 600)
    
    # Plot of fisher exact test
    output$fisher <- renderPlot({
        plot_results(compute_fisher(),input$plot_type)
    }, height = 600)
    
    # Plot of EWCE
    output$ewce <- renderPlot({
        plot_results(compute_EWCE(),input$plot_type)
        }, height = 600)
    
    # Plot heatmap
    output$heatmap <- renderPlot({
        plot_heatmap(dataset(),gene_list())
    }, height = 600)
    
    #Download table regression
    output$save_res_reg <- downloadHandler(
        filename = function() {
            datasetname <- parse_dataset_name(input$select)
            paste(input$file1,datasetname,"reg.txt",sep = ".")
        },
        content = function(file) {
            write_tsv(compute_reg(), file)
        }
    )
    
    #Download plot regression
    output$save_plot_reg <- downloadHandler(
        filename = function() {
            datasetname <- parse_dataset_name(input$select)
            paste(input$file1,datasetname,"reg.pdf",sep = ".")
        },
        content = function(file) {
            p <- plot_results(compute_reg(),input$plot_type)
            ggsave(file, plot = p, device = "pdf")
        }
    )
    
    #Download table fisher
    output$save_res_fisher <- downloadHandler(
        filename = function() {
            datasetname <- parse_dataset_name(input$select)
            paste(input$file1,datasetname,"fisher.txt",sep = ".")
        },
        content = function(file) {
            write_tsv(compute_fisher(), file)
        }
    )
    
    #Download plot fisher
    output$save_plot_fisher <- downloadHandler(
        filename = function() {
            datasetname <- parse_dataset_name(input$select)
            paste(input$file1,datasetname,"fisher.pdf",sep = ".")
        },
        content = function(file) {
            p <- plot_results(compute_fisher(),input$plot_type)
            ggsave(file, plot = p, device = "pdf")
        }
    )
    
    #Download table EWCE
    output$save_res_ewce <- downloadHandler(
        filename = function() {
            datasetname <- parse_dataset_name(input$select)
            paste(input$file1,datasetname,"ewce.txt",sep = ".")
        },
        content = function(file) {
            write_tsv(compute_EWCE(), file)
        }
    )
    
    #Download plot EWCE
    output$save_plot_ewce <- downloadHandler(
        filename = function() {
            datasetname <- parse_dataset_name(input$select)
            paste(input$file1,datasetname,"ewce.pdf",sep = ".")
        },
        content = function(file) {
            p <- plot_results(compute_EWCE(),input$plot_type)
            ggsave(file, plot = p, device = "pdf")
        }
    )
})
