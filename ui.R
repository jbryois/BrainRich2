# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.

# Author: Julien Bryois
# Date: 14.11.2019

library(shiny)
library(shinythemes)

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
            h5("An example gene set can be found",a("here",href="https://raw.githubusercontent.com/jbryois/BrainRich2/master/Additional_data/cahoy_oligodendrocyte.txt")),
            
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
                        choices = list("Zeisel et al. (2018) lvl2 (Mouse)" = "Data/Zeisel.lvl2.1to1.norm.txt.gz",
                                       "Zeisel et al. (2018) lvl3 (Mouse)" = "Data/Zeisel.lvl3.1to1.norm.txt.gz",
                                       "Zeisel et al. (2018) lvl4 (Mouse)" = "Data/Zeisel.lvl4.1to1.norm.txt.gz",
                                       "Zeisel et al. (2018) lvl5 (Mouse)" = "Data/Zeisel.1to1.norm.txt.gz",
                                       "Skene et al. (2018) lvl1 (Mouse)" = "Data/Skene_lvl1.1to1.norm.txt.gz",
                                       "Skene et al. (2018) lvl2 (Mouse)" = "Data/Skene_lvl2.1to1.norm.txt.gz",
                                       "Saunders et al. (2018) lvl1 (Mouse)" = "Data/Saunders.lvl1.1to1.norm.txt.gz",
                                       "Saunders et al. (2018) lvl2 (Mouse)" = "Data/Saunders.lvl2.1to1.norm.txt.gz",
                                       "Saunders et al. (2018) lvl3 (Mouse)" = "Data/Saunders.lvl3.1to1.norm.txt.gz",
                                       "Habib et al. (2017) (Human)" = "Data/Habib.norm.txt.gz",
                                       "GTex v7 (Human tissues)" = "Data/GTEx.v7.all.norm.txt.gz", 
                                       "GTex v8 (Human tissues)" = "Data/GTEx.v8.all.norm.txt.gz",
                                       "Allen Brain M1 (Human)" = "Data/AB_human_m1_10x.norm.txt.gz",
                                       "Allen Brain Multiple Cortical areas (Human)" = "Data/AB_multiple_cortical_areas_smartseq2019.norm.txt.gz",
                                       "Allen Brain MTG (Human)" = "Data/AB_mtg2018.norm.txt.gz",
                                       "Allen Brain Whole Cortex + hippocampus (Mouse)" = "Data/AB_whole_cortex_hippocampus_mouse_2020.norm.txt.gz"
                                       ),
                        selected = "Data/Zeisel.lvl4.1to1.norm.txt.gz"),
            
            # Input: Select plot type
            radioButtons("plot_type", label = h3("Plot Type"),
                         choices = list("-log10(pvalue)" = "-log10(pvalue)","Effect Sizes" = "Effect Sizes"), 
                         selected = "-log10(pvalue)")
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
                                 )
                        ),
                        tabPanel("Heatmap", fluid = TRUE,
                                 mainPanel(
                                     h5("Heatmap of the specificity metrics by decile."),
                                     h5("1 = lowest decile, 10 = top decile"),
                                     h5("X-axis: genes in the gene set"),
                                     plotOutput("heatmap")      
                                 )
                        ),
                        tabPanel("References", fluid = TRUE,
                                 mainPanel(
                                     h3("Datasets:"),
                                     h5(a("Zeisel et al. 2018",href="https://www.cell.com/cell/fulltext/S0092-8674(18)30789-X")),
                                     h5(a("Skene et al. 2018",href="https://www.nature.com/articles/s41588-018-0129-5")),
                                     h5(a("Habib et al. 2017",href="https://www.nature.com/articles/nmeth.4407")),
                                     h5(a("Saunders et al. 2018",href="https://www.sciencedirect.com/science/article/pii/S0092867418309553")),
                                     h5(a("GTEx v7",href="https://www.nature.com/articles/nature24277")),
                                     h5(a("GTEx v8",href="https://www.biorxiv.org/content/10.1101/787903v1")),
                                     h5(a("Allen Brain M1 (Human)",href="https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x")),
                                     h5(a("Allen Brain Multiple Cortical areas (Human)",href="https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq")),
                                     h5(a("Allen MTG (Human)",href="https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-smart-seq")),
                                     h5(a("Allen Brain Whole Cortex + hippocampus (Mouse)",href="https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-smart-seq")),
                                     h3("Methods:"),
                                     h5(a("EWCE",href="https://www.frontiersin.org/articles/10.3389/fnins.2016.00016/full")),
                                     h3("Code:"),
                                     h5(a("Github",href="https://github.com/jbryois/BrainRich2"))
                                 )
                        )
            )
        )
    )
)