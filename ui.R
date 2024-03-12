#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(DT)
#Load Data:

sgRNA_Tx_data<-read.table(file='./Data/GlycoCRISPR_sgRNALocs_ProteinCodingTx.txt',sep = '\t',header=FALSE)
colnames(sgRNA_Tx_data)<-c('GeneTx','sgRNA_id','sgRNA_start','sgRNA_end',
                           'TxLen','fracTarget','TxType','UniProtID')

#Get list of glycogenes to select:
lts_ggenes<-sapply(as.character(sgRNA_Tx_data$sgRNA_id),function(x) unlist(strsplit(x=x,split='\\_'))[4]) %>% 
  unique()

load('./Data/CCLE_cellTypes_Lines_List_NCI60Subset.rda')

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Glycogene CRISPR Library Viewer"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = 'geneSelect',label = 'Select Glycogene:',choices = lts_ggenes),
      uiOutput(outputId='Tx_Select_UI'),
      selectInput(inputId='tissue_select',label='Select Tissue Type:',choices=names(CCLE_cellline_types_NCI60)),
      uiOutput(outputId='cell_select'),
      actionButton('clearCells','Clear Cells'),
      textOutput('selected_cells'),br(),
      h3('Usage Instructions:'),
      h4('Downloads:'),
      tags$ul(
        tags$li(downloadLink('gcInstruct','Stepwise methods for using glycoCRISPR')),
        tags$li(downloadLink('mnscp','Author version of manuscript')),
        tags$li(downloadLink('suppTab','Supplementary Tables (includes gene targets and guide sequences)'))
      ),
      tags$span(style='color:red',"If you use this resource, please cite:\n"),
      
      tags$span(style='color:red','Y. Zhu, T. Groth, A. Kelkar, Y. Zhou, S. Neelamegham “A GlycoGene CRISPR-Cas9 lentiviral library to study lectin binding and human glycan biosynthesis pathways”, 
        Glycobiology, 2020 (in press) [.pdf].'),
      tags$em(''),
      br(),
      tags$span(style='color:blue',"All plasmids will be available from Addgene. If you need reagents urgently, please write to neel@buffalo.edu."),p(),
      
      p('This app allows users to browse sgRNA and associated data collected using the Glycogene CRISPR sgRNA library. To use:'),
      
      tags$ul(
        tags$li('"Select Glycogene:" using  drop-down menu'), 
        tags$li('“Select Gene Transcript:” as many genes have more than one protein-coding transcript variant.')
      ),
      
      p('Information about sgRNA targeting the selected genes, 
        location on transcript and experimental data regarding sgRNA enrichment/depletion will appear in table format.'),
      p('To determine glycogene transcripts appearing in selected, common cancer cell lines:'),
      
      tags$ul(
        tags$li('Select tissue type'), 
        tags$li('Select cell line')
      ),
      
      p('Data shown is rendered as heat maps using RSEM (RNA-seq by Expectation Maximization) collected from the Cancer Cell Line Encyclopedia (CCLE) (ref 3).
         An RSEM value of 1 corresponds to the top 20% of RSEM values in the CCLE for the selected NCI60 cell lines.  Thus, RSEM values > 1 
         are considered to be highly expressed.'),
      h4('References:'),
      tags$ol(
        tags$li('Y. Zhu, et al. Glycobiology, 2020'),
        tags$a(href="https://doi.org/10.1093/glycob/cwaa074", "[link]"),
        tags$li('Li, W. et al. Genome Biol 15, 554 (2014).'),
        tags$a(href="https://doi.org/10.1186/s13059-014-0554-4", "[link]"),
        tags$li('Barretina, J. et al. Nature 483, 603–607 (2012).'),
        tags$a(href="https://doi.org/10.1038/nature11003", "[link]")
        )
    ),
      
    mainPanel(
      h2('sgRNA Sequence and Location Table',align='center'),
      dataTableOutput('guide_infoTable'),
      h2('Gene Rank & False Discover Rate',align='center'),
      plotOutput("gene_RRA_Plot"),
      h2('sgRNA Log2 Fold Change Upon Lectin Enrichment'),
      dataTableOutput('guide_statTable'),
      h2('Selected CCLE Transcript RSEM Data',align='center'),
      plotOutput("gene_CCLEPlot")
    )
  )
))
