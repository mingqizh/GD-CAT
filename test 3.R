library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(reticulate)
library(Rcpp)
library(Seurat)
library(reshape2)
library(colormap)
library(bnlearn)
library(bnstruct)
library(patchwork)
library(dplyr)
library(limma)
library(MetBrewer)
library(qgraph)
library(devtools)
library(ADAPTS)
library(preprocessCore)
library(pheatmap)
library(ggplot2)
library(WGCNA)
library(mclust)
library(pheatmap)
library(qvalue)
library(DBI)
library(odbc)
library(RODBC)
library(RSQLite)
library(RPostgres)
library(enrichR)
ui<-fluidPage(
  tags$h1("Welcome to the GTEx Shiny App!"),
  br(),
  tags$img(height = 100, 
           width = 100, 
           src = "gtex.png"),
  sidebarLayout(
    sidebarPanel(
      radioButtons(
        inputId = "specie",
        label = "Chose a specie",
        choices = c("Human","Mouse")
      )
    ),
    mainPanel(verbatimTextOutput("specie"))
  ),
  br(),
  sidebarLayout(
    sidebarPanel(
      radioButtons(
        "Gender", "Chose Sex",c("Male", "Female", "Both")
      )
    ),
    mainPanel(verbatimTextOutput("gender"))
  ),
  br(),
  sidebarLayout(
    sidebarPanel(
      textInput("gene", "Input a gene"),
      textInput("origan_tissue", "Input a origan - tissue")
    ),
    mainPanel(
      verbatimTextOutput("gene"),
      verbatimTextOutput("origan_tissue")
    )),
  tags$h2("Let's see the pie"),
  br(),
  sidebarLayout(
    actionBttn(
      inputId = "bttn1",
      label = "Start",
      color = "primary",
      style = "bordered"
    ),
    plotOutput("res_bttn1_plot")
  ),
  plotOutput("res_bttn2_plot"),
  tableOutput("res_buttn1_table"),
  hr(),
  selectInput(
    "tissue2",
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
      Lung = "Lung", 
      Muscle_Skeletal = "Muscle - Skeletal", 
      Spleen = "Spleen", 
      Small_Intestine_Terminal_Ileum = "Small Intestine - Terminal Ileum", 
      Stomach = "Stomach", 
      Thyroid = "Thyroid", 
      Pancreas = "Pancreas", 
      Pituitary = "Pituitary")
  ),
  tabsetPanel(
    tabPanel("Pathways", 
             selectInput(
               "network",
               "Select the network",
               c(Undirected_network="Undirected network", 
                 Directed_network ="Directed network")
             ),
             plotOutput("plot3"),
             plotOutput("plot4"),
             plotOutput("plot5")),
    tabPanel("Enrichment", 
             tags$h4("GO_Biological_Process_2021"),
             plotOutput("enrichment1"),
             tags$h4("GO_Molecular_Function_2021"),
             plotOutput("enrichment2"),
             tags$h4("Reactome_2022"),
             plotOutput("enrichment3"),
             tags$h4("DSigDB"),
             plotOutput("enrichment4"))
  
))

server<-function(input, output){
  output$specie <- renderPrint({
    input$specie
  })
  output$gender <- renderPrint({
    input$Gender
  })
  working_data<-reactive({
    a<-odbcConnect("UCI_Seldin_lab", uid = "root", pwd = "1234567890zZ!")
    annot = sqlQuery(a, "USE gtex;", stringsAsFactors=F)
    annot = sqlQuery(a, "select * from a", stringsAsFactors=F)
    annot<-as.data.frame(annot)
  })
  output$gene <- renderPrint({
    input$gene
  })
  output$origan_tissue<-renderPrint({
    input$origan_tissue
  })

  gene_tissue<-reactive({
    origin_tissue = paste0(input$gene, input$origin_tissue)
  })
 

  output$res_bttn1_plot = renderPlot({
    annot<-working_data()
    sig_table = annot[annot$qvalue<0.1,]
    sig_table$qcat =ifelse(sig_table$qvalue<0.01, 'q<0.01', 'q<0.1')
    sig_table$qcat =ifelse(sig_table$qvalue<0.0001, 'q<0.0001', paste0(sig_table$qcat))
    table(sig_table$tissue[sig_table$qcat=='q<0.01'])
    
    blank_theme <- theme_minimal()+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=25, face="bold")
      )
    
    binned_sig_prots= sig_table %>%
      dplyr::group_by(qcat, tissue_2) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::mutate(freq = n / sum(n))
    
    
    tissue_freqs = binned_sig_prots %>% dplyr::group_by(qcat) %>% dplyr::summarise(sum(n))
    
    binned_sig_prots$tot_count = tissue_freqs$`sum(n)`[match(binned_sig_prots$qcat, tissue_freqs$qcat)]
    binned_sig_prots$qcat1 = paste0(binned_sig_prots$qcat, ' ', binned_sig_prots$tot_count, ' genes')
    col_scheme = rev(met.brewer('Austria', length(unique(sig_table$tissue_2))))
    names(col_scheme) = unique(sig_table$tissue)
    binned_sig_prots$tissue_2 <-  reorder(binned_sig_prots$tissue_2, binned_sig_prots$n)
    

    ggplot(binned_sig_prots, aes(x = "", y = freq, fill =tissue_2)) + 
      geom_bar(stat = "identity", width = 1, position = position_fill()) + blank_theme + theme(plot.title=element_text(size=25, face="bold")) +
      theme(axis.text.x=element_blank())+ scale_fill_manual(values=col_scheme) +
      coord_polar(theta = "y") + 
      guides(fill = guide_legend(reverse = TRUE)) +
      facet_wrap( ~ qcat1)


  })
  output$res_bttn2_plot = renderPlot({
    annot<-working_data()
    annot<-annot[!grepl(annot$tissue_1[1], annot$tissue_2),]
    sig_table = annot[annot$qvalue<0.1,]
    sig_table$qcat =ifelse(sig_table$qvalue<0.01, 'q<0.01', 'q<0.1')
    sig_table$qcat =ifelse(sig_table$qvalue<0.0001, 'q<0.0001', paste0(sig_table$qcat))
    table(sig_table$tissue[sig_table$qcat=='q<0.01'])
    
    blank_theme <- theme_minimal()+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=25, face="bold")
      )
    
    binned_sig_prots= sig_table %>%
      dplyr::group_by(qcat, tissue_2) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::mutate(freq = n / sum(n))
    
    
    tissue_freqs = binned_sig_prots %>% dplyr::group_by(qcat) %>% dplyr::summarise(sum(n))
    
    binned_sig_prots$tot_count = tissue_freqs$`sum(n)`[match(binned_sig_prots$qcat, tissue_freqs$qcat)]
    binned_sig_prots$qcat1 = paste0(binned_sig_prots$qcat, ' ', binned_sig_prots$tot_count, ' genes')
    col_scheme = rev(met.brewer('Austria', length(unique(sig_table$tissue_2))))
    names(col_scheme) = unique(sig_table$tissue)
    binned_sig_prots$tissue_2 <-  reorder(binned_sig_prots$tissue_2, binned_sig_prots$n)
    
    
    ggplot(binned_sig_prots, aes(x = "", y = freq, fill =tissue_2)) + 
      geom_bar(stat = "identity", width = 1, position = position_fill()) + blank_theme + theme(plot.title=element_text(size=25, face="bold")) +
      theme(axis.text.x=element_blank())+ scale_fill_manual(values=col_scheme) +
      coord_polar(theta = "y") + 
      guides(fill = guide_legend(reverse = TRUE)) +
      facet_wrap( ~ qcat1)
    
    
  })
  # pathway network
  all_tog2<-reactive({
    new_table <- working_data()
    
    new_table = new_table[new_table$pvalue<1e-4,]
    new_table = new_table[!is.na(new_table$tissue),]
    sig_set = new_table[!new_table$tissue==gene_tissue,]
    row_length = ifelse(length(row.names(sig_set))>300, as.numeric(300), as.numeric(paste0(length(row.names(sig_set)))))
    sig_set = sig_set[1:row_length,]
    
    all_tog2 = working_dataset[, colnames(working_dataset) %in% sig_set$gene_tissue | colnames(working_dataset)== gene_tissue]
    all_tog2[is.na(all_tog2)] = 0
  })
  new_cols<-reactive({
    working_dataset<-working_data()
    tissue_set = unique(gsub(".*_","",working_dataset$gene_tissue_2))
    new_cols = met.brewer('Moreau', length(tissue_set))
    names(new_cols) = unique(tissue_set)
  })
  colkey1<-reactive({
    all_tog2<-all_tog2()
    colkey1 = as.data.frame(gsub(".*_","",colnames(all_tog2)))
    colnames(colkey1) = 'tissue'
    
    colkey1$cols = new_cols[match(colkey1$tissue, names(new_cols))]
  })
  map2<-reactive({
    all_tog2<-all_tog2()
    bics_map = bicorAndPvalue(all_tog2, all_tog2, use = 'p')
    
    map1 = as.data.frame(bics_map$bicor)
    map1 = melt(as.matrix(map1))
    map1$value[map1$value > 0.999999] <- 0
    map2 = dcast(map1, Var1 ~ Var2, value.var = 'value')
    row.names(map2) = map2$Var1
    map2$Var1 = NULL
  })
  
  
  # plot undirected network
  output$plot3 = renderPlot({
    map2<-map2()
    colkey1<-colkey1()
    p<-qgraph(map2, minimum = 0.2, cut = 0.6, vsize = 3, color=colkey1$cols, legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=2, directed=F, labels = colnames(map2)) + ggtitle('')
    print(p)
  })
  
  #plot directed network
  output$plot4 = renderPlot({
    all_tog2<-all_tog2()
    colkey1<-colkey1()
    net1 = mmhc(all_tog2) 
    p<-qgraph(net1, vsize = 3, legend = F, color=colkey1$cols, borders = TRUE, layout='spring', label.cex=2, directed=T, labels = colnames(map2)) + ggtitle('')
    p<-as.data.frame(p)
    print(p)
  })
  
  # plot tissue legend
  output$plot5 = renderPlot({
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    legend("topleft", legend =names(new_cols), pch=16, pt.cex=3, cex=1.5, bty='n',
           col = new_cols)
    mtext("Tissue", at=0.2, cex=2)
  })
  
  # enrichment 
  enriched = reactive({
    annot <- working_data()

    sig_table = annot[annot$qvalue<0.1,]
    sig_table$qcat =ifelse(sig_table$qvalue<0.01, 'q<0.01', 'q<0.1')
    sig_table$qcat =ifelse(sig_table$qvalue<0.0001, 'q<0.0001', paste0(sig_table$qcat))
    table(sig_table$qcat)
    
    binned_sig_prots= sig_table %>%
      dplyr::group_by(tissue_2, qcat) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::mutate(freq = n / sum(n))
    
    
    tissue_list = binned_sig_prots[binned_sig_prots$qcat=='q<0.01',]
    tissue_list = tissue_list[order(tissue_list$n, decreasing = T),]
    
    select_tissue = input$tissue2
    pp1 = annot[annot$qvalue<0.05,]
    pp1 = pp1[pp1$tissue_2 %in% select_tissue,]
    pp1_length = ifelse(length(row.names(pp1)) > 200, as.numeric(200), as.numeric(length(row.names(pp1))))
    pp2 = pp1[1:pp1_length,]
    gg1 = pp2$gene_symbol_2
    
    setEnrichrSite("Enrichr")
    dbs <- listEnrichrDbs()
    dbs1 <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "Reactome_2022", "DSigDB")
    
    enriched <- enrichr(gg1, dbs1)
    

  })

  output$enrichment1<-renderPlot({
    enriched<-enriched()
    plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
  })
  output$enrichment2<-renderPlot({
    enriched<-enriched()
    plotEnrich(enriched[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
  })
  output$enrichment3<-renderPlot({
    enriched<-enriched()
    plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
  })
  output$enrichment4<-renderPlot({
    enriched<-enriched()
    plotEnrich(enriched[[4]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
  })
}

shinyApp(ui = ui, server = server)

ui <- dashboardPage(
  +     dashboardHeader(),
  +     dashboardSidebar(),
  +     dashboardBody(plotOutput("res_bttn1_plot")))


a<-odbcConnect("UCI_Seldin_lab", uid = "root", pwd = "1234567890zZ!")
annot = sqlQuery(a, "USE gtex;", stringsAsFactors=F)
annot = sqlQuery(a, "select * from first100gene", stringsAsFactors=F)
