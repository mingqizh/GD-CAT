library(shiny)
library(shinydashboard)
library(dplyr)
library(purrr)
library(ggplot2)
library(enrichR)
library(echarts4r)
library(MetBrewer)
library(forcats)
library(DT)
library(zip)

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
          fluidPage(
            box(title="Top-N",status = "primary",width=12,collapsible = TRUE,solidHeader = TRUE,
                sidebarLayout(
                  sidebarPanel(width=3,
                               radioButtons(inputId = "specie",label = "Chose a specie",choices = c("Human","Mouse")),
                               radioButtons("Gender", "Chose Sex",c("Male", "Female", "Both")),
                               sliderInput(inputId = "topn",label = "How many genes that correlated with the tissue do you want to see?",value = 30, min = 1, max = 50)),
                  mainPanel(plotOutput('plot.top'))
                )),
            tabBox(
              title = "Pie Chart",
              width = 12,
              tabPanel('Pie-1', echarts4rOutput('plot.p1')),
              tabPanel('Pie-2', echarts4rOutput('plot.p2')),
              tabPanel('Pie-3', echarts4rOutput('plot.p3')),
              tabPanel('Table-1',DTOutput('table.t1')),
              tabPanel('Pie-4', echarts4rOutput('plot.p4')),
              tabPanel('Pie-5', echarts4rOutput('plot.p5')),
              tabPanel('Pie-6', echarts4rOutput('plot.p6')),
              tabPanel('Table-2',DTOutput('table.t2'))
            ),
            box(
              title = "Pie Chart",
              width = 12,
              verbatimTextOutput('text'),
              actionButton('btn', class = "btn-primary", 'Enrichment Analysis'),
              uiOutput('plots.en')
            )
          )
  )
}

body <- function(){
  dashboardBody(
    tabItems(
      t2()
    )
  )
}

working_data<- function(){
  read.csv('D:/R-lab/shiny-pie-click/test.csv')
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

  data<-binned_sig_prots
  echart<- binned_sig_prots %>%
    group_by(qcat) %>%
    group_split() %>%
    map(~{
      qcat1<-unique(.x$qcat1)
      qcat<-sub('q<','',unique(.x$qcat))
      .x %>%
        e_chart(tissue_2) %>%
        e_pie(freq,name = qcat,right='20%') %>%
        e_title(qcat1,left="center",textStyle=list(fontSize=12)) %>%
        e_legend(type ='scroll',
                 orient='vertical',
                 top='center',
                 right='5%') %>%
        e_on(
          list(seriesName = qcat),
          "function(x){
          //alert(Object.keys(x));
          var msg = [x.seriesIndex,x.seriesName,x.name,x.dataIndex]
          //alert(msg)
          Shiny.setInputValue('selected_tissue',x.name, {priority: 'event'})
          Shiny.setInputValue('selected_q',x.seriesName, {priority: 'event'});
        }"
        ) %>%
        e_toolbox_feature(feature = "saveAsImage",title='Save')
    })
  return(list(data=data,echart=echart))
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
  data<-binned_sig_prots # data
  echart<-binned_sig_prots %>% # chart
    group_by(qcat) %>%
    group_split() %>%
    map(~{
      qcat1<-unique(.x$qcat1)
      qcat<-sub('q<','',unique(.x$qcat))
      .x %>%
        e_chart(tissue_2) %>%
        e_pie(freq,name = qcat,right='20%') %>%
        e_title(qcat1,left="center",textStyle=list(fontSize=12)) %>%
        e_legend(type ='scroll',
                 orient='vertical',
                 top='center',
                 right='5%') %>%
        e_on(
          list(seriesName = qcat),
          "function(x){
          //alert(Object.keys(x));
          var msg = [x.seriesIndex,x.seriesName,x.name,x.dataIndex]
          //alert(msg)
          Shiny.setInputValue('selected_tissue',x.name, {priority: 'event'})
          Shiny.setInputValue('selected_q',x.seriesName, {priority: 'event'});
        }"
        ) %>%
        e_toolbox_feature(feature = "saveAsImage",title='保存')
    })
  return(list(data=data,echart=echart))
}

f_topgen <- function(working_data,max){
  col_scheme= rev(met.brewer('Austria', length(unique(working_data$tissue_2))))
  names(col_scheme) = unique(working_data$tissue_2)
  max_gene_length = max
  max_gene_length<-as.numeric(max_gene_length)
  top_genes = working_data[1:max_gene_length,]
  top_genes$color = col_scheme[match(top_genes$tissue_2, names(col_scheme))]

  ggplot(top_genes, aes(x=fct_reorder(tissue_2, abs(bicor), .desc=T), y=bicor, fill=tissue_2)) + geom_col(position = position_dodge2(reverse = TRUE)) +  theme_classic() + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values = top_genes$color[order(abs(top_genes$bicor), decreasing = T)]) + xlab('') + ylab('bicor coefficent') + ggtitle('Top genes correlated with origin gene tissue')
}

f_e1 <- function(annot,select_tissue,q){
  pp1 = annot[annot$qvalue<as.numeric(q),]
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

f_e2 <- function(annot,select_tissue,q){
  pp1 = annot[annot$qvalue<as.numeric(q),]
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
  # pie chart
  pie1<-f1(working_data())$echart
  output$plot.p1 <- renderEcharts4r({
    pie1[[1]]
  })
  output$plot.p2 <- renderEcharts4r({
    pie1[[2]]
  })
  output$plot.p3 <- renderEcharts4r({
    pie1[[3]]
  })
  pie2<-f2(working_data())$echart
  output$plot.p4 <- renderEcharts4r({
    pie2[[1]]
  })
  output$plot.p5 <- renderEcharts4r({
    pie2[[2]]
  })
  output$plot.p6 <- renderEcharts4r({
    pie2[[3]]
  })
  # table
  output$table.t1 <- renderDT(
        f1(working_data())$data,
        extensions = 'Buttons',
        options = list(dom = 'Blfrtip',
                       scrollX=TRUE,
                       autoWidth=FALSE,
                       buttons = c('copy','csv','excel'),
                       lengthMenu = list(c(10.25,50),
                                        c(10,25,59,"All"))))
  output$table.t2 <- renderDT(
    f2(working_data())$data,
    extensions = 'Buttons',
    options = list(dom = 'Blfrtip',
                   scrollX=TRUE,
                   autoWidth=FALSE,
                   buttons = c('copy','csv','excel'),
                   lengthMenu = list(c(10.25,50),
                                     c(10,25,59,"All"))
    ))
  # top n
  output$plot.top<- renderPlot(
    f_topgen(working_data(),input$topn)
  )
  # click
  output$text <- renderText({
    print(c(input$selected_tissue,input$selected_q))
  })
  # enrichment
  plots<-reactiveVal()
  observeEvent( input$btn, {
    if(!is.null(input$selected_tissue)){
      progress <- Progress$new(session, min=0, max=5)
      on.exit(progress$close())
      progress$set(message = 'Calculation in progress',
                   detail = 'This may take a while...')
      progress$set(value = 1)

      enriched1<-f_e1(working_data(),input$selected_tissue,input$selected_q)
      plots_en1<- enriched1 %>%
        map(~plotEnrich(.x, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value"))
      progress$set(value = 2)

      enriched2<-f_e2(working_data(),input$selected_tissue,input$selected_q)
      plots_en2<- enriched2 %>%
        map(~plotEnrich(.x, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value"))
      progress$set(value = 3)

      plots_all<-c(plots_en1,plots_en2)
      plots(plots_all)
      output[['plots.en']]<-renderUI({
        plot_output_list <- lapply(1:length(plots_all), function(i) {
          plotname <- paste("en", i, sep="")
          plotOutput(plotname)
        })
        plot_output_list$btn_down<- downloadButton("download", "Download Image")
        do.call(tagList, plot_output_list)
      }
      )
      progress$set(value = 4)
      for (i in 1:length(plots_all)) {
        local({
          my_i <- i
          plotname <- paste("en", my_i, sep="")
          output[[plotname]] <- renderPlot({
            plots_all[[my_i]]
          })
        })
      }
      progress$set(value = 5)
    } else {
      showNotification("Please first click the Pie chart body to selected an tissue.",duration = NULL,type="warning")
    }
  })
  output$download <- downloadHandler(
    filename = function() {
      paste("Plots-", Sys.Date(), ".zip", sep="")
    },
    content = function(file) {
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      plots<-plots()
      for(i in 1:length(plots)){
        ggsave( paste0('plot',i,'.png'), plot = plots[[i]], device = "png")
      }
      zip::zip(file,paste0('plot',1:length(plots),'.png'))
    }
  )
  # tips
  observeEvent(input$tabs, {
    if(input$tabs=='t2'){
      showNotification("Click the legend (on the right of the chart body) of the Pie chart to toggle the display of the series.",duration = NULL,type="message")
    }
  })
}

shinyApp(ui, server)


