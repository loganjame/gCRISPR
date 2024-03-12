#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(dplyr)
library(stringr)
library(gridExtra)
library(DT)

# ggplot(test %>% filter(id=='ST6GAL1')) + geom_hline(yintercept = 0,size=2,color='black') + geom_point(aes(x=pos.rank,y=y,colour=group),size=3) + geom_label_repel(aes(x=pos.rank,y=y,label=apply(test %>% filter(id=='ST6GAL1'),1,function(x) paste(x['group'],paste('Enrichment FDR =',x['pos.fdr']),paste('Rank =',x['pos.rank']),sep = '\n')),fill=group),force=6) + labs(y='',x='Gene Ranking by RRA') + theme_light() +theme(panel.border = element_blank(),panel.grid=element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = 'none')  
#Load Data:
#load('./Data/sgRNA_Tx_data_SriramLFC.rda')
sgRNA_Tx_data_LFC<-read.table(file='./Data/sgRNA_Tx_data_SriramLFC.tsv',sep='\t',header=TRUE)
sgRNA_Tx_data_LFC<-sgRNA_Tx_data_LFC %>% mutate(
  TxId=sapply(as.character(GeneTx),function(x) unlist(strsplit(x=x,split='\\|'))[2]),
  fracTarget=sapply(fracTarget,function(x) round(x,3))
)
#CCLE Transcript data:

# CCLE_data<-read.table(file='./Data/CCLE_RNAseq_rsem_transcripts_tpm_NCI60Subset.txt',sep='\t',header=TRUE)
#load('./Data/CCLE_glycoData.rda')
CCLE_glycoData<-read.table(file='./Data/CCLE_glycoData.tsv',sep='\t',header=TRUE)
CCLE_glycoData$transcript_id<-sapply(as.character(CCLE_glycoData$transcript_id),function(x){
  unlist(strsplit(split='\\.',x=x))[1]
})
#Cell tissue --> types
load('./Data/CCLE_cellTypes_Lines_List_NCI60Subset.rda')

#Old RRA data:
#load('./Data/Total_RRA_CRISPRdata.rda')
#load('./Data/Total_RRA_CRISPR_sigGroups.rda')

#RRA analysis:
load('./Data/totalNorm_results_merged.rda')
totalNorm_results_merged=as.data.frame(totalNorm_results_merged)
totalNorm_results_merged$group=as.character(totalNorm_results_merged$group)
#Color Ramp:
cols<-colorRampPalette(c('blue','white','red'))

colorVec<-c(cols(1000))
brks<-seq(-10,10,20/1000)[-1:-2]

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  output$gcInstruct<-downloadHandler(filename='GlycogeneCRISPR_instructions.pdf',
    content=function(con){
      file.copy('GlycogeneCRISPR_instructions.pdf',con)
    }
  )
  output$suppTab<-downloadHandler(filename='SupplementalTables_060620_sub.xls',
    content=function(con){
      file.copy('SupplementalTables_060620_sub.xls',con)
    }
  )
  output$mnscp<-downloadHandler(filename='Manuscript_AuthorVersion.pdf',
    content=function(con){
      file.copy('Manuscript_AuthorVersion.pdf',con)
    }
  )
  
  #Process cell selection logic:
  output$cell_select<-renderUI({
    selectInput('cell_select','Select Cell Type:',choices=c("",CCLE_cellline_types_NCI60[[input$tissue_select]]),
                selected = NULL)
  })
  #Process Tx selection info:
  output$Tx_Select_UI<-renderUI({
    Tx_List<-sgRNA_Parse<-sgRNA_Tx_data_LFC %>% filter(sapply(sgRNA_id,function(x) grepl(input$geneSelect,x))==TRUE) %>%
      .[,'TxId'] %>% unique()
    selectInput('Tx_selection','Select Gene Transcript:',choices = c("",Tx_List))
  })
  
  cellList<-reactiveVal(c())
  
  observeEvent(input$cell_select,{
    newList<-c(cellList(),input$cell_select)
    cellList(unique(newList))
  })
  
  observeEvent(input$clearCells,{
    cellList(c())
  })
  
  output$selected_cells<-renderText({
    paste('Selected Cells:',
          paste(cellList(),collapse = '\n'),sep = '\n')
  })
  #End cell selection logic
  
  output$guide_infoTable <- renderDataTable({
    #Parse the gene selection:
    sgRNA_Parse<-sgRNA_Tx_data_LFC %>% filter(Gene==input$geneSelect) %>% 
      filter(TxId==input$Tx_selection) %>%
      dplyr::select(-TxId,-Gene,-TxLen,-fracTarget,-TxType) %>% dplyr::select(-contains('results'))
    dtaDT<-datatable(sgRNA_Parse %>% 
            rename(
              `sgRNA ID`=sgRNA_id,
              `sgRNA Sequence`=sgRNA_seq,
              `Gene and Transcript`=GeneTx,
              `sgRNA Transcript Start`=sgRNA_start,
              `sgRNA Transcript End`=sgRNA_end,
              `UniProt ID`=UniProtID)
            ,options = list(scrollX=TRUE,dom='<p<t>l>',pageLength=min(c(15,dim(sgRNA_Parse)[1])),rownames = FALSE))
    
    
  })
  
  output$gene_RRA_Plot<-renderPlot({
    #New Total normalization-based RRA data below:
    RRAParse<-totalNorm_results_merged %>% filter(id==input$geneSelect) %>% arrange(pos.rank)
    RRAParse$y<-rep(0,dim(RRAParse)[1])
    # New Total Normalization-based data has new "group" names, created a new
    # mapping:
    mapList=list(
      'E_minus'='E-Selectin Low',
      'P_minus'='P-Selectin Low',
      'HECA_minus'='HECA-452 Low',
      'PNA_plus'='Peanut agglutinin (PNA) High',
      'PHAL_minus'='Phaseolus Vulgaris Leucoagglutinin (PHA-L) Low',
      'VVA_plus'='Vicia Villosa Lectin (VVA) High'
    )
    RRAParse$group<-sapply(RRAParse$group,function(x) mapList[[x]])
    message(RRAParse)
    plt<-ggplot(RRAParse) + geom_hline(yintercept = 0,size=2,color='black') + 
      geom_label_repel(aes(x=pos.rank,y=y,label=apply(RRAParse,1,
        function(x) paste(x['group'],paste('Enrichment FDR =',round(as.numeric(x['pos.fdr']),3)),paste('Rank =',x['pos.rank']),sep = '\n')),fill=group),
        force=75,direction = 'both',min.segment.length = 150,size=5) + geom_point(aes(x=pos.rank,y=y),size=7,color='black') + 
      geom_point(aes(x=pos.rank,y=y,colour=group),size=6) + 
      scale_x_continuous(breaks=seq(0,length(unique(RRAParse$id)),50)) + expand_limits(x=c(0,350)) + 
      labs(y='',x='Gene Ranking by Robust Rank Aggregation (ref 2)') + theme_light() + 
      theme(panel.border = element_blank(),panel.grid=element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),
            legend.position = 'none',axis.title.x=element_text(size=12))
    return(plt)
  },height=350)
  
  output$guide_statTable <- renderDataTable({
    #Parse the gene selection:
    sgRNA_Parse<-sgRNA_Tx_data_LFC %>% filter(Gene==input$geneSelect) %>% 
      filter(TxId==input$Tx_selection) %>% distinct() %>% dplyr::select(sgRNA_id,contains('results')) %>% 
      rename(`E-Selectin Low`=E_minus_results,
             `P-Selectin Low`=P_minus_results,
             `HECA-452 Low`=HECA_minus_results,
             `Peanut agglutinin (PNA) High`=PNA_plus_results,
             `Phaseolus Vulgaris Leucoagglutinin (PHA-L) Low`=PHAL_minus_results,
             `Vicia Villosa Lectin (VVA) High`=VVA_plus_results)
    # sgRNA_Parse<- sgRNA_Parse %>% mutate_at(names(sgRNA_Parse) %>% .[grepl('results',.)],function(x) factor(x,levels=c('---','--','-','+','++','+++','NS')))
    
    dtaDT<-datatable(sgRNA_Parse,options = list(scrollX=TRUE,dom='<p<t>l>',pageLength=min(c(15,dim(sgRNA_Parse)[1]))),rownames = FALSE) %>%
      formatStyle(names(sgRNA_Parse) %>% .[grepl('(High|Low)',.)],backgroundColor = styleInterval(brks,colorVec))
    return(dtaDT)
    
  })

  output$gene_CCLEPlot <- renderPlot({
    #Parse the gene selection:
    sgRNA_Parse<-sgRNA_Tx_data_LFC %>% filter(sapply(sgRNA_id,function(x) grepl(input$geneSelect,x))==TRUE)
    
    if (length(cellList())==0){
      CCLE_parse<-NULL
    } else {
      CCLE_parse<-CCLE_glycoData %>% filter(transcript_id %in% sgRNA_Parse$TxId) %>%
        dplyr::select(gene_id,transcript_id,one_of(cellList()))
    }
    
    
    #Create expression grid:
    if (is.null(CCLE_parse)){
      return(NULL)
    }
    if (dim(CCLE_parse)[2]==2){
      return(NULL)
    } else {
      CCLE_parse_melt<-melt(CCLE_parse,id.vars=c('gene_id','transcript_id'))
      CCLE_parse_melt$variable<-sapply(as.character(CCLE_parse_melt$variable),function(x) str_wrap(gsub('\\_','\\ ',x=x),25))
      p_CCLE_grid<-ggplot(CCLE_parse_melt) + geom_tile(aes(x=transcript_id,y=variable,fill=value),color='black',size=1.75) + 
        geom_text(aes(x=transcript_id,y=variable,label=signif(value,3)),size=6) + 
        labs(x='Transcript',y='Cell Line',title=paste(input$geneSelect,'Transcript Abundances')) +
        scale_fill_gradientn(colours = c('white','red','red'),breaks=c(0,0.5,1)) + 
        coord_flip() + theme_light() + theme(axis.text.x=element_text(size=12,hjust=1,angle=45),axis.text.y=element_text(size=12),
                                             axis.title=element_text(size=20),title=element_text(size=15),
                                             panel.border = element_blank(),panel.grid = element_blank())
    }
    
    #Return the CCLE grid:
    return(p_CCLE_grid)
    
  })
  
})
