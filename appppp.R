library(shiny)
library(shinydashboard)
library(shinyalert)
library(dplyr)
library(purrr)
library(ggplot2)
library(enrichR)
library(echarts4r)
library(MetBrewer)
library(forcats)
library(DT)
library(zip)
library(reticulate)
library(Rcpp)
library(reshape2)
library(colormap)
library(patchwork)
library(qgraph)
library(devtools)
library(preprocessCore)
library(pheatmap)
library(WGCNA)
library(mclust)
library(pheatmap)
library(qvalue)
library(RPostgres)
library(enrichR)
library(stats)
library(forcats)
library(MetBrewer)
library(feather)
library(tidyr)
library(tibble)
library(webshot)
allowWGCNAThreads()

header <- function(){
  dashboardHeader(title = "Welcome to GD-CAT!")
}

sidebar <- function(){
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      class= "hide",
      menuItem("Instruction",  tabName = "t2",icon = icon("dashboard")),
      menuItem("Settings",  tabName = "t2",icon = icon("dashboard"))
    ),
    radioButtons(inputId = "specie",label = "Chose a species",choices = c("Human","Mouse - coming soon")),
    radioButtons("Gender", "Chose Sex",c("Both"="both","Male"="male", "Female"="female")),
    textInput("origin_gene","origin_gene, Official NBCI gene symbol",value = "ADIPOQ"),
    selectInput(
      "origin_tissue",
      "Select a tissue",
      c(Adipose_Subcutaneous="Adipose - Subcutaneous",
        Adipose_Visceral = "Adipose - Visceral (Omentum)", 
        Adrenal_Gland = "Adrenal Gland", 
        Artery_Aorta = "Artery - Aorta", 
        Artery_Coronary = "Artery - Coronary", 
        Brain_Hippocampus = "Brain - Hippocampus", 
        Brain_Hypothalamus = "Brain - Hypothalamus", 
        Colon_Transverse = "Colon - Transverse", 
        Colon_Sigmoid = "Colon - Sigmoid", 
        Heart_Left_Ventricle = "Heart - Left Ventricle", 
        Kidney_Cortex = "Kidney - Cortex", 
        Liver = "Liver", 
        Muscle_Skeletal = "Muscle - Skeletal", 
        Spleen = "Spleen", 
        Small_Intestine_Terminal_Ileum = "Small Intestine - Terminal Ileum", 
        Stomach = "Stomach", 
        Thyroid = "Thyroid", 
        Pancreas = "Pancreas", 
        Pituitary = "Pituitary")
    ),
    actionButton('import', class = "btn-primary", 'Process data')
  )
}

t2<-function(){
  tabItem(tabName="t2",
          fluidPage(
            tabBox(
              title = "Pie Chart",
              width = 12,
              tabPanel('All genes q<0.1', echarts4rOutput('plot.p1')),
              tabPanel('All genes q<0.01', echarts4rOutput('plot.p2')),
              tabPanel('All genes q<0.001', echarts4rOutput('plot.p3')),
              tabPanel('Origin-removed q<0.1', echarts4rOutput('plot.p4')),
              tabPanel('Origin-removed q<0.01', echarts4rOutput('plot.p5')),
              tabPanel('Origin-removed q<0.001', echarts4rOutput('plot.p6')),
              tabPanel('Table-1',DTOutput('table.t1'))
            ),
            shinydashboard::box(
              title = "Enrichement Analysis",
              width = 12,
              verbatimTextOutput('text'),
              actionButton('btn', class = "btn-primary", 'Start Analysis'),
              uiOutput('plots.en')
            ),
            tabBox(title="Cell type",width=12,height="500px",
                   tabPanel('Cell type',echarts4rOutput('cell'))
            ),
            sliderInput(inputId = "topn",label = "How many top-ranked correlated genes do you want to see?",value = 30, min = 1, max = 50),
            tabBox(title="Top-N",width=12,height="500px",
                   tabPanel('Top-ranked genes',echarts4rOutput('plot.top1')),
                   tabPanel('Top-ranked genes without origin tissue',echarts4rOutput('plot.top2'))
            ),
            sliderInput(inputId = "within",label = "Within tissue gene numbers",value = 30, min = 1, max = 50),
            sliderInput(inputId = "external",label = "Peripheral tissue gene numbers",value = 100, min = 1, max = 150),
            shinydashboard::box(
              title = "Network Analysis",
              width = 12,
              actionButton('btn1', class = "btn-primary", 'Start Analysis'),
              tabPanel('Chart',plotOutput('net')),
              downloadButton("PDFPlot", "Download PDF"),
              tabPanel('Chart',echarts4rOutput('plot.map3')),
              tabPanel('Table',DTOutput('table.map2'))
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

ui<-dashboardPage(header(), sidebar(), body())

f1<-function(annot, col_scheme){
  
  
  binned_sig_prots= annot %>%
    dplyr::group_by(qcat, tissue_2) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(freq = n / sum(n))%>%
    dplyr::arrange(desc(freq))
  
  
  tissue_freqs = binned_sig_prots %>% dplyr::group_by(qcat) %>% dplyr::summarise(sum(n))
  
  binned_sig_prots$tot_count = tissue_freqs$`sum(n)`[match(binned_sig_prots$qcat, tissue_freqs$qcat)]
  binned_sig_prots$qcat1 = paste0(binned_sig_prots$qcat, ' ', binned_sig_prots$tot_count, ' genes')
  
  data<-binned_sig_prots
  binned_sig_prots$color = col_scheme[match(binned_sig_prots$tissue_2, names(col_scheme))]
  echart<- binned_sig_prots %>%
    group_by(qcat) %>%
    group_split() %>%
    purrr::map(~{
      qcat1<-unique(.x$qcat1)
      qcat<-sub('q<','',unique(.x$qcat))
      .x %>%
        e_chart(tissue_2) %>%
        e_pie(freq,name = qcat,right='20%') %>%
        e_add_nested("itemStyle", color) %>%
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

f2<-function(annot,origin_tissue, col_scheme){
  annot<-annot[!grepl(origin_tissue, annot$tissue_2),]
  sig_table = annot[annot$qvalue<0.1,]
  sig_table1 = annot[annot$qvalue<0.01,]
  sig_table1$qcat = paste0(0.01)
  sig_table2 = annot[annot$qvalue<0.001,]
  sig_table2$qcat = paste0(0.001)
  sig_table = as.data.frame(rbind(sig_table, sig_table1, sig_table2))
  
  binned_sig_prots= sig_table %>%
    dplyr::group_by(qcat, tissue_2) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(freq = n / sum(n))%>%
    dplyr::arrange(desc(freq))
  
  
  tissue_freqs = binned_sig_prots %>% dplyr::group_by(qcat) %>% dplyr::summarise(sum(n))
  
  binned_sig_prots$tot_count = tissue_freqs$`sum(n)`[match(binned_sig_prots$qcat, tissue_freqs$qcat)]
  binned_sig_prots$qcat1 = paste0(binned_sig_prots$qcat, ' ', binned_sig_prots$tot_count, ' genes')
  binned_sig_prots$tissue_2 <-  reorder(binned_sig_prots$tissue_2, binned_sig_prots$n)
  data<-binned_sig_prots # data
  binned_sig_prots$color = col_scheme[match(binned_sig_prots$tissue_2, names(col_scheme))]
  echart<-binned_sig_prots %>% # chart
    group_by(qcat) %>%
    group_split() %>%
    purrr::map(~{
      qcat1<-unique(.x$qcat1)
      qcat<-sub('q<','',unique(.x$qcat))
      .x %>%
        e_chart(tissue_2) %>%
        e_pie(freq,name = qcat,right='20%') %>%
        e_add_nested("itemStyle", color) %>%
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
        e_toolbox_feature(feature = "saveAsImage",title='ä¿å­')
    })
  return(list(data=data,echart=echart))
}

f_e1 <- function(annot,select_tissue,q){
  pp1 = annot[annot$qvalue<as.numeric(q),]
  pp1 = pp1[pp1$bicor>0,]
  pp1 = pp1[pp1$tissue_2 %in% select_tissue,]
  pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
  pp2 = pp1[1:pp1_length,]
  gg1 = pp2$gene_symbol_2
  
  setEnrichrSite("Enrichr")
  dbs <- listEnrichrDbs()
  dbs1 <- c("GO_Biological_Process_2021", "Reactome_2022")
  
  enriched <- enrichr(gg1, dbs1)
  
}

f_e2 <- function(annot,select_tissue,q){
  pp1 = annot[annot$qvalue<as.numeric(q),]
  pp1 = pp1[pp1$bicor<0,]
  pp1 = pp1[pp1$tissue_2 %in% select_tissue,]
  pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
  pp2 = pp1[1:pp1_length,]
  gg1 = pp2$gene_symbol_2
  
  setEnrichrSite("Enrichr")
  dbs <- listEnrichrDbs()
  dbs1 <- c("GO_Biological_Process_2021", "Reactome_2022")
  
  enriched <- enrichr(gg1, dbs1)
  
  
}


get_cell<-function(working_dataset, sig_table, origin_gene, origin_tissue,col_scheme){
  #read in sc-seq matrix 
  sc_matrix = read.csv('C:/Users/mingqiz7/Desktop/GTEx app/data/full pan-tissue DECON abundances.csv')
  origin_gene_tissue = paste0(origin_gene, '_', origin_tissue)
  #filter sc-seq matrix based on working_dataset object 
  row.names(sc_matrix) = sc_matrix$GTExID
  sc_matrix$GTExID = NULL
  sc_matrix = sc_matrix[row.names(sc_matrix) %in% row.names(working_dataset),]
  tissue1 <- working_dataset[,colnames(working_dataset) %in% origin_gene_tissue]
  
  full_cors_sc = bicorAndPvalue(tissue1, sc_matrix, use = 'p')
  cor_table_sc = reshape2::melt(full_cors_sc$bicor)
  new_p = reshape2::melt(full_cors_sc$p)
  cor_table_sc$Var1=NULL
  colnames(cor_table_sc) = c('cell_tissue', 'bicor')
  
  
  cor_table_sc$pvalue = signif(new_p$value, 3)
  cor_table_sc$bicor = round(cor_table_sc$bicor, 3)
  cor_table_sc$qvalue = signif(p.adjust(cor_table_sc$pvalue, "BH"), 3)
  cor_table_sc = cor_table_sc[order(cor_table_sc$qvalue, decreasing=F),]
  cor_table_sc = na.omit(cor_table_sc)
  cor_table_sc$tissue = gsub(".*_","",cor_table_sc$cell_tissue)
  cor_table_sc$tissue = gsub('...', ' - ', cor_table_sc$tissue, fixed=T)
  cor_table_sc$tissue = gsub('..Omentum.', ' (Omentum)', cor_table_sc$tissue, fixed=T)
  cor_table_sc$tissue = gsub('Left.Ventricle', 'Left Ventricle', cor_table_sc$tissue, fixed=T)
  top_genes2 = cor_table_sc
  top_genes2$color = col_scheme[match(top_genes2$tissue, names(col_scheme))]
  top_genes2 = top_genes2[order(abs(top_genes2$bicor), decreasing = T),]
}
get_top_genes1<-function(sig_table,max_gene_length, col_scheme){
  top_genes = sig_table[1:max_gene_length,]
  top_genes$color = col_scheme[match(top_genes$tissue_2, names(col_scheme))]
  top_genes
}
get_top_genes2<-function(sig_table,max_gene_length,origin_tissue, col_scheme){
  
  sig_table1 = sig_table[!sig_table$tissue_2 %in% origin_tissue,]
  top_genes = sig_table1[1:max_gene_length,]
  top_genes$color = col_scheme[match(top_genes$tissue_2, names(col_scheme))]
  
  top_genes
}

server <- function(input, output, session) {
  working_dataset<-eventReactive(input$import,{
    load("C:/Users/mingqiz7/Desktop/GTEx app/data/working_dataset.RData")
    isolate({
      if(input$Gender == "both"){
        working_dataset<-both
      }else if(input$Gender == "male"){
        working_dataset<-male
      }else{
        working_dataset<-female
      }
      working_dataset
    })
  })

  
  sig_table<-eventReactive(input$import,{
    progress <- Progress$new(session, min=0, max=5)
    on.exit(progress$close())
    progress$set(message = 'Pre-processing raw data',
                 detail = 'It will take around 2-5 minutes depending on usage')
    progress$set(value = 1)
    isolate({
      origin_gene = input$origin_gene
      origin_tissue = input$origin_tissue
      origin_gene_tissue = paste0(origin_gene, '_', origin_tissue)
      working_dataset<-working_dataset()
    })
    
    tissue2 <- working_dataset[,grepl('Adipose - Subcutaneous', colnames(working_dataset)) | grepl('Adipose - Visceral (Omentum)', colnames(working_dataset), fixed=T) | grepl('Brain - Hypothalamus', colnames(working_dataset)) | grepl('Brain - Hippocampus', colnames(working_dataset)) | grepl('Small Intestine - Terminal Ileum', colnames(working_dataset), fixed=T) | grepl('Stomach', colnames(working_dataset), fixed=T) | grepl('Thyroid', colnames(working_dataset), fixed=T) | grepl('Pancreas', colnames(working_dataset), fixed=T) | grepl('Spleen', colnames(working_dataset), fixed=T) | grepl('Muscle - Skeletal', colnames(working_dataset), fixed=T) | grepl('Pituitary', colnames(working_dataset), fixed=T) | grepl('Artery - Coronary', colnames(working_dataset), fixed=T) | grepl('Liver', colnames(working_dataset), fixed=T) | grepl('Kidney - Cortex', colnames(working_dataset), fixed=T) | grepl('Heart - Left Ventricle', colnames(working_dataset), fixed=T) | grepl('Colon - Transverse', colnames(working_dataset), fixed=T) | grepl('Colon - Sigmoid', colnames(working_dataset), fixed=T) | grepl('Adrenal Gland', colnames(working_dataset), fixed=T) |  grepl('Artery - Aorta', colnames(working_dataset), fixed=T),]
    
    
    tissue1 <- working_dataset[,colnames(working_dataset) %in% origin_gene_tissue]
    
    
    #tissue1 = tissue1[row.names(tissue1) %in% row.names(tissue2),]
    #tissue2 = tissue2[row.names(tissue2) %in% row.names(tissue1),]
    progress$set(value = 2)
    full_cors = bicorAndPvalue(tissue1, tissue2, use = 'p')
    progress$set(value = 3)
    cor_table = reshape2::melt(full_cors$bicor)
    new_p = reshape2::melt(full_cors$p)
    
    colnames(cor_table) = c('gene_tissue_1', 'gene_tissue_2', 'bicor')
    #can drop here to clear CPU
    full_cors=NULL
    
    cor_table$pvalue = signif(new_p$value, 3)
    cor_table$bicor = round(cor_table$bicor, 3)
    cor_table = cor_table[!cor_table$gene_tissue_1==cor_table$gene_tissue_2,]
    cor_table$qvalue = signif(p.adjust(cor_table$pvalue, "BH"), 3)
    cor_table = cor_table[order(cor_table$qvalue, decreasing=F),]
    cor_table = na.omit(cor_table)
    cor_table$gene_symbol_1 = gsub("\\_.*","",cor_table$gene_tissue_1)
    cor_table$tissue_1 = gsub(".*_","",cor_table$gene_tissue_1)
    cor_table = cor_table[!is.na(cor_table$tissue_1),]
    progress$set(value = 4)
    cor_table$gene_symbol_2 = gsub("\\_.*","",cor_table$gene_tissue_2)
    cor_table$tissue_2 = gsub(".*_","",cor_table$gene_tissue_2)
    cor_table = cor_table[!is.na(cor_table$tissue_2),]
    
    #sig_table åé¢éè¦ç¨å°
    sig_table = cor_table[cor_table$qvalue<0.1,]
    sig_table$qcat =ifelse(sig_table$qvalue<0.01, 'q<0.01', 'q<0.1')
    sig_table$qcat =ifelse(sig_table$qvalue<0.001, 'q<0.001', paste0(sig_table$qcat))
    
    sig_table
    
  })
  tryCatch({
    col_scheme<-reactive({
    sig_table<-sig_table()
    col_scheme = rev(met.brewer('Signac', length(unique(sig_table$tissue_2))))
    names(col_scheme) = unique(sig_table$tissue_2)
    col_scheme
  }) 
  }, warning = function(w) {
    shinyalert("Warning!","Message to explain this", type = "warning")
  })
  
  #cell type
  output$cell<- renderEcharts4r({
    top_genes1<-get_cell(working_dataset(),sig_table(), input$origin_gene, input$origin_tissue, col_scheme())
    
    top_genes1 %>%
      arrange(desc(abs(bicor))) %>%
      e_chart(cell_tissue, width='400px', height = NULL) %>%
      e_bar(bicor) %>%
      e_legend(show=F) %>%
      e_add_nested("itemStyle", color) %>%
      e_grid(bottom="150px") %>%
      e_x_axis(axisLabel = list(interval = 0, rotate = 45)) %>%
      e_y_axis(name = "Bicor Values")%>%
      e_datazoom(type='inside') %>%
      e_toolbox(show=F) %>%
      e_tooltip(
        trigger = 'item',
        axisPointer = list(
          type = "shadow",
          axis='x'
        )
      )
  }
  )
  # top n
  output$plot.top1<- renderEcharts4r({
    top_genes1<-get_top_genes1(sig_table(),input$topn, col_scheme())
    
    top_genes1 %>%
      arrange(desc(abs(bicor))) %>%
      e_chart(gene_tissue_2, width='400px', height = NULL) %>%
      e_bar(bicor) %>%
      e_legend(show=F) %>%
      e_add_nested("itemStyle", color) %>%
      e_grid(bottom="150px") %>%
      e_x_axis(axisLabel = list(interval = 0, rotate = 45)) %>%
      e_y_axis(name = "Bicor Values")%>%
      e_datazoom(type='inside') %>%
      e_toolbox(show=F) %>%
      e_tooltip(
        trigger = 'item',
        axisPointer = list(
          type = "shadow",
          axis='x'
        )
      )
  }
  )
  output$plot.top2<- renderEcharts4r({
    top_genes2<-get_top_genes2(sig_table(),input$topn,input$origin_tissue, col_scheme())
    
    top_genes2 %>%
      arrange(desc(abs(bicor))) %>%
      e_chart(gene_tissue_2, width='500px', height = NULL) %>%
      e_bar(bicor) %>%
      e_legend(show=F) %>%
      e_add_nested("itemStyle", color) %>%
      e_grid(bottom="100px") %>%
      e_x_axis(axisLabel = list(interval = 0, rotate = 45)) %>%
      e_y_axis(name = "Bicor Values")%>%
      e_datazoom(type='inside') %>%
      e_toolbox(show=F) %>%
      e_tooltip(
        trigger = 'item',
        axisPointer = list(
          type = "shadow",
          axis='x'
        )
      )
  }
  )
 
  net<-eventReactive(input$btn1,{
    progress <- Progress$new(session, min=0, max=5)
    on.exit(progress$close())
    progress$set(message = 'Generating the enrichements',
                 detail = 'Just wait a second')
    progress$set(value = 1)
    isolate({
      number_orig_gene =as.numeric(input$within) 
    number_peripheral_genes =as.numeric(input$external) 
    sig_table<-sig_table()
    working_dataset<-working_dataset()
    origin_gene<-input$origin_gene
    origin_tissue=input$origin_tissue
    })
    
    origin_gene_tissue = paste0(origin_gene, '_', origin_tissue)
    col_scheme<- col_scheme()
    origin_pull = sig_table[sig_table$tissue_2==origin_tissue,]
    orig_network_genes = as.vector(origin_pull$gene_tissue_2[1:number_orig_gene])
    peripheral_pull = sig_table[!sig_table$tissue_2==origin_tissue,]
    periph_network_genes = as.vector(peripheral_pull$gene_tissue_2[1:number_peripheral_genes])
    sql_pull_list = c(origin_gene_tissue, orig_network_genes, periph_network_genes)
    
    progress$set(value = 2)
    
    tissue1 <- working_dataset[, colnames(working_dataset)  %in% sql_pull_list]
    full_cors = bicorAndPvalue(tissue1, tissue1, use = 'p')
    cor_table = reshape2::melt(full_cors$bicor)
    new_p = reshape2::melt(full_cors$p)
    colnames(cor_table) = c('gene_tissue_1', 'gene_tissue_2', 'bicor')
    
    progress$set(value = 3)
    #can drop here to clear CPU
    full_cors=NULL
    cor_table$pvalue = signif(new_p$value, 3)
    cor_table$bicor = round(cor_table$bicor, 3)
    cor_table = cor_table[!cor_table$gene_tissue_1==cor_table$gene_tissue_2,]
    cor_table$qvalue = signif(p.adjust(cor_table$pvalue, "BH"), 3)
    cor_table = cor_table[order(cor_table$qvalue, decreasing=F),]
    cor_table = na.omit(cor_table)
    cor_table$gene_symbol_1 = gsub("\\_.*","",cor_table$gene_tissue_1)
    cor_table$tissue_1 = gsub(".*_","",cor_table$gene_tissue_1)
    cor_table = cor_table[!is.na(cor_table$tissue_1),]
    cor_table$gene_symbol_2 = gsub("\\_.*","",cor_table$gene_tissue_2)
    cor_table$tissue_2 = gsub(".*_","",cor_table$gene_tissue_2)
    cor_table = cor_table[!is.na(cor_table$tissue_2),]
    sql_pull_2 = cor_table
    
    progress$set(value = 4)
    ######################################
    #the new SQL pull table is not listed as sql_pull_2
    #this will take work on your end to interface with pulling the correct genes from tissue_1 based on the number selected in the app.  I arbitrarily chose 10 here and you obviously wouldn't jsut apply top_n since it will differ.  The final list is called network_gene_lsit
    map1 = sql_pull_2
    map1$cols = col_scheme[match(map1$tissue_1, names(col_scheme))]
    map2 = dcast(map1, gene_tissue_1 ~ gene_tissue_2, value.var = 'bicor', fun.aggregate = mean)
    row.names(map2) = map2$gene_tissue_1
    map2$gene_tissue_1 = NULL
    cols_set = as.data.frame(row.names(map2))
    row.names(cols_set) = row.names(map2)
    cols_set$cols = map1$cols[match(row.names(cols_set), map1$gene_tissue_1)]
    cols_set$cols = ifelse(row.names(cols_set) %in% origin_gene_tissue, 'gray5', paste0(cols_set$cols))
    
    qgraph(map2, minimum = 0.5, cut = 0.85, vsize = 3, color=cols_set$cols, legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=4, directed=F, labels = colnames(map2)) 
  })
  output$net<-renderPlot({
    net()
  })
  output$PDFPlot <- downloadHandler(
    filename = function() {
      paste("Net-", Sys.Date(), ".pdf", sep="")
    },
    content = function(file) {
        pdf(file)
        plot(net())
        dev.off()
    }
  )
  #output$table.map2<-renderDT(
  #  get_map2(working_dataset(),sig_table(),input$slicen,input$origin_gene_tissue),
  # extensions = 'Buttons',
  # options = list(dom = 'Blfrtip',
  #                scrollX=TRUE,
  #                autoWidth=FALSE,
  #                buttons = c('copy','csv','excel'),
  #                lengthMenu = list(c(10.25,50),
  #                                  c(10,25,59,"All"))
  # )
  #)
  # pie chart
  pie1<-reactive({
    print(dim(sig_table()))
    f1(sig_table(), col_scheme())$echart
  })
  output$plot.p1 <- renderEcharts4r({
    pie1()[[3]]
  })
  output$p1 <- downloadHandler(
    filename = function() {
      paste("Pie-", Sys.Date(), ".pdf", sep="")
    },
    content = function(file) {
      webshot(pie1()[[3]]$html, file)
    }
  )
  output$plot.p2 <- renderEcharts4r({
    pie1()[[2]]
  })
  output$plot.p3 <- renderEcharts4r({
    pie1()[[1]]
  })
  pie2<-reactive({
    f2(sig_table(),input$origin_tissue, col_scheme())$echart
  })
  output$plot.p4 <- renderEcharts4r({
    pie2()[[2]]
  })
  output$plot.p5 <- renderEcharts4r({
    pie2()[[1]]
  })
  output$plot.p6 <- renderEcharts4r({
    pie2()[[3]]
  })
  
  # table
  table1<-reactive({
    sig_table<-sig_table()
    sig_table<-sig_table[,-c(1, 6, 7)]
  })
  output$table.t1 <- renderDT(
    table1(),
    extensions = 'Buttons',
    options = list(dom = 'Blfrtip',
                   scrollX=TRUE,
                   autoWidth=FALSE,
                   buttons = c('copy','csv','excel'),
                   lengthMenu = list(c(10.25,50),
                                     c(10,25,59,"All"))))
  
  
  # click
  output$text <- renderText({
    print(c('You selected tissue:',input$selected_tissue,', and q<',input$selected_q))
  })
  # enrichment
  plots<-reactiveVal()
  observeEvent( input$btn, {
    if(!is.null(input$selected_tissue)){
      progress <- Progress$new(session, min=0, max=5)
      on.exit(progress$close())
      progress$set(message = 'Gnerating the enrichement',
                   detail = 'Just wait a second')
      progress$set(value = 1)
      tryCatch({
        # Your plot code here
        enriched1<-f_e1(sig_table(),input$selected_tissue,input$selected_q)
      plots_en1<- enriched1 %>%
        purrr::map(~plotEnrich(.x, showTerms = 10, numChar = 30, y = "Count", orderBy = "P.value") 
                   + ggtitle(paste0('Positive gene correlations with ', input$origin_gene, ' ',input$origin_tissue, ' ', names(.x))))
      progress$set(value = 2)
      
      enriched2<-f_e2(sig_table(),input$selected_tissue,input$selected_q)
      plots_en2<- enriched2 %>%
        purrr::map(~plotEnrich(.x, showTerms = 10, numChar = 30, y = "Count", orderBy = "P.value")
                   + ggtitle(paste0('Negative gene correlations with ', input$origin_gene, ' ', input$origin_tissue, ' ',names(.x))))
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
      }, error = function(e) {
        # Display an error message if an error occurs
        shinyalert("Oops","The selceted tissue do not contain enough gene to generate the enrichment. Please select again.", type = "error")
       })
      
      progress$set(value = 5)
    } else {
      shinyalert("Warning!", "Please first click the Pie chart body to selected an tissue.", type = "warning")
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
      showNotification("Click the tissue on the Pie chart to toggle the display of the series.",duration = NULL,type="message")
    }
  })
}

shinyApp(ui, server)


