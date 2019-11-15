# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.

# Author: Julien Bryois
# Date: 14.11.2019

library(shiny)
library(shinythemes)
library(tidyverse)
library(readxl)
library(broom)

# Define UI for application
ui <- fluidPage(
    
    # App title
    titlePanel("BrainRich2"),
    
    # Theme 
    theme = shinytheme("flatly"),                   
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            h4(strong("Gene Set Enrichment in Brain Cell Types")),
            h5("Please input a file with your gene list (first column)."),
            h5("The app works with ensembl, symbol or entrez gene ids."),
            h5("An example gene set can be found",a("here",href="https://raw.githubusercontent.com/jbryois/BrainRich2/master/03_Example/cahoy_oligodendrocyte.txt")),
            
            # Input: Select type of delimiter for file to be uploaded
            radioButtons('sep', 'Separator',
                         c(Tab='\t',
                           Space=' ',
                           Comma=',',
                           xls='xls',
                           xlsx='xlsx'),
                         '\t'),
            # Input: Check whether the file has a header or not
            checkboxInput('header', 'Header', TRUE),
            
            # Input: Load input file
            fileInput('file1', 'Input File'),
            
            # Input: Select Dataset to use
            selectInput("select", label = h3("Select Dataset"), 
                        choices = list("Zeisel et al. (2018) lvl2" = "Data/Zeisel.lvl2.1to1.norm.txt.gz",
                                       "Zeisel et al. (2018) lvl3" = "Data/Zeisel.lvl3.1to1.norm.txt.gz",
                                       "Zeisel et al. (2018) lvl4" = "Data/Zeisel.lvl4.1to1.norm.txt.gz",
                                       "Zeisel et al. (2018) lvl5" = "Data/Zeisel.1to1.norm.txt.gz",
                                       "Skene et al. (2018) lvl1" = "Data/Skene_lvl1.1to1.norm.txt.gz",
                                       "Skene et al. (2018) lvl2" = "Data/Skene_lvl2.1to1.norm.txt.gz",
                                       "Saunders et al. (2018) lvl1" = "Data/Saunders.lvl1.1to1.norm.txt.gz",
                                       "Saunders et al. (2018) lvl2" = "Data/Saunders.lvl2.1to1.norm.txt.gz",
                                       "Saunders et al. (2018) lvl3" = "Data/Saunders.lvl3.1to1.norm.txt.gz",
                                       "Habib et al. (2017)" = "Data/Habib.norm.txt.gz",
                                       "GTex v7" = "Data/GTEx.v7.all.norm.txt.gz", 
                                       "GTex v8" = "Data/GTEx.v8.all.norm.txt.gz"),
                        selected = "Data/Zeisel.lvl4.1to1.norm.txt.gz"),
            
            # Input: Select plot type
            radioButtons("plot_type", label = h3("Plot Type"),
                         choices = list("-log10(pvalue)" = "-log10(pvalue)","Effect Sizes" = "Effect Sizes"), 
                         selected = "-log10(pvalue)"),
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
            
            # Output: Tabset w/ plot ----
            tabsetPanel(type = "tabs",
                        tabPanel("Enrichment Regression", fluid = TRUE,
                                mainPanel(
                                    h5("Tests whether the average specificity of the gene list is greater than 0 using linear regression."),
                                    h5("The specificity metrics are z-scaled for each tissue/cell type prior to regression"),
                                    downloadButton('save_res_reg',label = "Save Table"),
                                    downloadButton('save_plot_reg',label = "Save Plot"),
                                    plotOutput("reg")      
                                    )
                        ),
                                 
                        tabPanel("Enrichment Fisher", fluid = TRUE,
                                 mainPanel(
                                     h5("Tests whether the gene list is enriched among the 10% most specific genes in each tissue/cell type."),
                                     downloadButton('save_res_fisher',label = "Save Table"),
                                     downloadButton('save_plot_fisher',label = "Save Plot"),
                                     plotOutput("fisher")      
                                 )
                        ),
                        tabPanel("Enrichment EWCE", fluid = TRUE,
                                 mainPanel(
                                     h5("Tests whether the average specificity of the gene list is greater than expected by chance"),
                                     h5("Selects 10'000 random gene lists of the same size to estimate enrichment."),
                                     strong("Please be patient, this may take a while."),
                                     br(),
                                     downloadButton('save_res_ewce',label = "Save Table"),
                                     downloadButton('save_plot_ewce',label = "Save Plot"),
                                     plotOutput("ewce")      
                                 )),
                        tabPanel("Heatmap", fluid = TRUE,
                                 mainPanel(
                                     h5("Heatmap of the specificity metrics by decile."),
                                     h5("1 = lowest decile, 10 = top decile"),
                                     h5("X-axis: genes in the gene set"),
                                     plotOutput("heatmap")      
                                 )),
                        tabPanel("References", fluid = TRUE,
                                 mainPanel(
                                     h3("Datasets:"),
                                     h5(a("Zeisel et al. 2018",href="https://www.cell.com/cell/fulltext/S0092-8674(18)30789-X")),
                                     h5(a("Skene et al. 2018",href="https://www.nature.com/articles/s41588-018-0129-5")),
                                     h5(a("Habib et al. 2017",href="https://www.nature.com/articles/nmeth.4407")),
                                     h5(a("Saunders et al. 2018",href="https://www.sciencedirect.com/science/article/pii/S0092867418309553")),
                                     h5(a("GTEx v7",href="https://www.nature.com/articles/nature24277")),
                                     h5(a("GTEx v8",href="https://www.biorxiv.org/content/10.1101/787903v1")),
                                     h3("Methods:"),
                                     h5(a("EWCE",href="https://www.frontiersin.org/articles/10.3389/fnins.2016.00016/full")),
                                     h3("Code:"),
                                     h5(a("Github",href="https://github.com/jbryois/BrainRich2"))
                                 )
                        )
            ),
        )
    )
)