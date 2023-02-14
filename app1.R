library(shiny)
library(shinydashboard)
library(dplyr)
library(purrr)
library(ggplot2)
library(enrichR)
library(echarts4r)
library(MetBrewer)
library(forcats)

header <- function(){
  dashboardHeader(title = "Welcome to the GTEx Shiny App!")
}

sidebar <- function(){
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Instruction", tabName = "t1", icon = icon("dashboard")),
      menuItem("Settings",  tabName = "t2",icon = icon("fa-solid fa-microscope"))
    )
  )
}

t2<-function(){
  tabItem(tabName="t2",
          fluidRow(
            column(width=3,
                   box(title="Predictiors Input",status = "primary",width=12,collapsible = TRUE,solidHeader = TRUE,
                       radioButtons(
                         inputId = "specie",
                         label = "Chose a specie",
                         choices = c("Human","Mouse")
                       ),
                       radioButtons(
                         "Gender", "Chose Sex",c("Male", "Female", "Both")
                       ),
                       textInput("gene", "Input a gene"),
                       textInput("origan_tissue", "Input a origan - tissue"),
                       sliderInput(inputId = "topn",
                                   label = "How many genes that correlated with the tissue do you want to see?",
                                   value = 30, min = 1, max = 50
                       ),
                       plotOutput('plot.top1'),
                       plotOutput('plot.top2'))),
            column(width=9,
                   tabBox(title="Pie Chart",width=12,
                          tabPanel('Pie-1',echarts4rOutput('plot.1')),
                          tabPanel('Pie-2',echarts4rOutput('plot.2'),
                                   verbatimTextOutput('clicked2')),
                          tabPanel('Pie-3',echarts4rOutput('plot.3'),
                                   verbatimTextOutput('clicked3'))),
                   box(title="Pie Chart",width=12,
                       verbatimTextOutput('text'),
                       actionButton('btn','Analysis'),
                       plotOutput('plot.e1'),
                       plotOutput('plot.e2'),
                       plotOutput('plot.e3')
                   ))
          ))
}

body <- function(){
  dashboardBody(
    tabItems(
      t2()
    )
  )
}

working_data<- function(){
  a<-odbcConnect("UCI_Seldin_lab", uid = "root", pwd = "1234567890zZ!")
  annot = sqlQuery(a, "USE gtex;", stringsAsFactors=F)
  annot = sqlQuery(a, "select * from a", stringsAsFactors=F)
  annot<-as.data.frame(annot)
}

f1<-function(annot){
  
  sig_table = annot[annot$qvalue<0.1,]
  sig_table$qcat =ifelse(sig_table$qvalue<0.01, 'q<0.01', 'q<0.1')
  sig_table$qcat =ifelse(sig_table$qvalue<0.0001, 'q<0.0001', paste0(sig_table$qcat))
  table(sig_table$tissue[sig_table$qcat=='q<0.01'])
  
  binned_sig_prots= sig_table %>%
    dplyr::group_by(qcat, tissue_2) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(freq = n / sum(n))%>%
    dplyr::arrange(desc(freq))
  
  
  tissue_freqs = binned_sig_prots %>% dplyr::group_by(qcat) %>% dplyr::summarise(sum(n))
  
  binned_sig_prots$tot_count = tissue_freqs$`sum(n)`[match(binned_sig_prots$qcat, tissue_freqs$qcat)]
  binned_sig_prots$qcat1 = paste0(binned_sig_prots$qcat, ' ', binned_sig_prots$tot_count, ' genes')
  
  binned_sig_prots %>%
    group_by(qcat1) %>%
    group_split() %>%
    map(~{
      qcat1<-unique(.x$qcat1)
      .x %>%
        e_chart(tissue_2) %>%
        e_pie(freq,name = qcat1,right='20%') %>%
        e_title(qcat1,left="center",textStyle=list(fontSize=12)) %>%
        e_legend(type ='scroll',
                 orient='vertical',
                 top='center',
                 left='80%') %>%
        e_on(
          list(seriesName = qcat1),
          "function(x){
          //alert(Object.keys(x));
          var msg = [x.seriesIndex,x.seriesName,x.name,x.dataIndex]
          alert(msg)
          Shiny.setInputValue('selected_tissue',x.name, {priority: 'event'});
        }"
        )
    })
  
}

f2<-function(annot){
  annot<-annot[!grepl(annot$tissue_1[1], annot$tissue_2),]
  sig_table = annot[annot$qvalue<0.1,]
  sig_table$qcat =ifelse(sig_table$qvalue<0.01, 'q<0.01', 'q<0.1')
  sig_table$qcat =ifelse(sig_table$qvalue<0.0001, 'q<0.0001', paste0(sig_table$qcat))
  
  binned_sig_prots= sig_table %>%
    dplyr::group_by(qcat, tissue_2) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(freq = n / sum(n))
  
  
  tissue_freqs = binned_sig_prots %>% dplyr::group_by(qcat) %>% dplyr::summarise(sum(n))
  
  binned_sig_prots$tot_count = tissue_freqs$`sum(n)`[match(binned_sig_prots$qcat, tissue_freqs$qcat)]
  binned_sig_prots$qcat1 = paste0(binned_sig_prots$qcat, ' ', binned_sig_prots$tot_count, ' genes')
  binned_sig_prots$tissue_2 <-  reorder(binned_sig_prots$tissue_2, binned_sig_prots$n)
  
  binned_sig_prots %>%
    group_by(qcat1) %>%
    group_split() %>%
    map(~{
      qcat1<-unique(.x$qcat1)
      .x %>%
        e_chart(tissue_2) %>%
        e_pie(freq,name = qcat1,right='20%') %>%
        e_title(qcat1,left="center",textStyle=list(fontSize=12)) %>%
        e_legend(type ='scroll',
                 orient='vertical',
                 top='center',
                 left='80%') %>%
        e_on(
          list(seriesName = qcat1),
          "function(x){
          alert(Object.keys(x));
          Shiny.setInputValue('chart_on_click',x.dataIndex +'-'+ x.name, {priority: 'event'});
        }"
        )
    })
}

f_topgen1 <- function(working_data,max){
  col_scheme= rev(met.brewer('Austria', length(unique(working_data$tissue_2))))
  names(col_scheme) = unique(working_data$tissue_2)
  max_gene_length = max
  max_gene_length<-as.numeric(max_gene_length)
  top_genes = working_data[1:max_gene_length,]
  top_genes$color = col_scheme[match(top_genes$tissue_2, names(col_scheme))]
  
  ggplot(top_genes, aes(x=fct_reorder(tissue_2, abs(bicor), .desc=T), y=bicor, fill=tissue_2)) + geom_col(position = position_dodge2(reverse = TRUE)) +  theme_classic() + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values = top_genes$color[order(abs(top_genes$bicor), decreasing = T)]) + xlab('') + ylab('bicor coefficent') + ggtitle('Top genes correlated with origin gene tissue')
}

f_topgen2 <- function(working_data,max){
  col_scheme= rev(met.brewer('Austria', length(unique(working_data$tissue_2))))
  names(col_scheme) = unique(working_data$tissue_2)
  max_gene_length = max
  max_gene_length<-as.numeric(max_gene_length)
  working_data = working_data[!working_data$tissue_2 %in% origin_tissue,]
  top_genes = working_data[1:max_gene_length,]
  top_genes$color = col_scheme[match(top_genes$tissue_2, names(col_scheme))]
  
  ggplot(top_genes, aes(x=fct_reorder(tissue_2, abs(bicor), .desc=T), y=bicor, fill=tissue_2)) + geom_col(position = position_dodge2(reverse = TRUE)) +  theme_classic() + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values = top_genes$color[order(abs(top_genes$bicor), decreasing = T)]) + xlab('') + ylab('bicor coefficent') + ggtitle('Top genes correlated with origin gene tissue')
}

f_e1 <- function(annot,select_tissue){
  pp1 = annot[annot$qvalue<0.05,]
  pp1 = pp1[pp1$bicor>0,]
  pp1 = pp1[pp1$tissue_2 %in% select_tissue,]
  pp1_length = ifelse(length(row.names(pp1)) > 200, as.numeric(200), as.numeric(length(row.names(pp1))))
  pp2 = pp1[1:pp1_length,]
  gg1 = pp2$gene_symbol_2
  
  setEnrichrSite("Enrichr")
  dbs <- listEnrichrDbs()
  dbs1 <- c("GO_Biological_Process_2021", "KEGG_2021_Human")
  
  enriched <- enrichr(gg1, dbs1)
  
}

f_e2 <- function(annot,select_tissue){
  pp1 = annot[annot$qvalue<0.05,]
  pp1 = pp1[pp1$bicor<0,]
  pp1 = pp1[pp1$tissue_2 %in% select_tissue,]
  pp1_length = ifelse(length(row.names(pp1)) > 200, as.numeric(200), as.numeric(length(row.names(pp1))))
  pp2 = pp1[1:pp1_length,]
  gg1 = pp2$gene_symbol_2
  
  setEnrichrSite("Enrichr")
  dbs <- listEnrichrDbs()
  dbs1 <- c("GO_Biological_Process_2021", "KEGG_2021_Human")
  
  enriched <- enrichr(gg1, dbs1)
  
  
}

#enriched1<-f_e1(working_data(),'Adipose - Subcutaneous')
#plotEnrich(enriched1[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

ui<-dashboardPage(header(), sidebar(), body())

server <- function(input, output, session) {
  pie1<-f1(working_data())
  output$plot.1 <- renderEcharts4r({
    pie1[[1]]
  })
  output$plot.2 <- renderEcharts4r({
    pie1[[2]]
  })
  output$plot.3 <- renderEcharts4r({
    pie1[[3]]
  })
  # top n
  output$plot.top1<- renderPlot(
    f_topgen1(working_data(),input$topn)
  )
  output$plot.top2<- renderPlot(
    f_topgen2(working_data(),input$topn)
  )
  # click
  output$text <- renderText({
    print(input$selected_tissue)
  })
  # enrichment
  observeEvent( input$btn, {
    if(!is.null(input$selected_tissue)){
      progress <- Progress$new(session, min=1, max=7)
      on.exit(progress$close())
      
      progress$set(message = 'Calculation in progress',
                   detail = 'This may take a while...')
      
      enriched1<-f_e1(working_data(),input$selected_tissue)
      progress$set(value = 1)
      
      output$plot.e1 <- renderPlot({
        plotEnrich(enriched1[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
      })
      progress$set(value = 2)
      output$plot.e2 <- renderPlot({
        plotEnrich(enriched1[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
      })
      progress$set(value = 3)
      output$plot.e3 <- renderPlot({
        plotEnrich(enriched1[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
      })
      progress$set(value = 4)
      enriched2<-f_e2(working_data(),input$selected_tissue)
      progress$set(value = 5)
      output$plot.e4 <- renderPlot({
        plotEnrich(enriched2[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
      })
      progress$set(value = 6)
      output$plot.e5 <- renderPlot({
        plotEnrich(enriched2[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
      })
      progress$set(value = 7)
      output$plot.e6 <- renderPlot({
        plotEnrich(enriched2[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
      })
    }
  })
  
}

shinyApp(ui, server)


