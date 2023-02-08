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
  library(highcharter)
  library(billboarder)
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
  if (interactive()) {

  header <- dashboardHeader(title = "Welcome to the GTEx Shiny App!")
  
  sidebar <- dashboardSidebar(
    sidebarUserPanel("Let's start it"),
    sidebarMenu(
      # Setting id makes input$tabs give the tabName of currently-selected tab
      id = "tabs",
      menuItem("Instruction", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Settings",  tabName = "Settings",icon = icon("fa-solid fa-microscope")),
      menuItem("Let's see the pie chart", tabName = "Pie",icon = icon("fa-doutone fa-chart-pie")),
      menuItem("The network", tabName = "pathway",icon = icon("fa-doutone fa-diagram-project")),
      menuItem("The enrichemtn",tabName = "enrich", icon = icon("fa-doutone fa-chart-bar"))
    )
  )
  
  body <- dashboardBody(
    tabItems(
      tabItem("dashboard",
              div(p("Dashboard tab content"))
      ),
      tabItem("Settings",
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
              br(),
              sidebarLayout(
                sidebarPanel(
                  textInput(inputId = "i", 
                              label = "How many genes that correlated with the tissue do you want to see?"
                              )
                ),
                mainPanel(
                  verbatimTextOutput('i'),
                  plotOutput('topgene1')
                ))
      ),
      tabItem("Pie",
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
              tableOutput("res_buttn1_table")
      ),
      tabItem("pathway",
              "Sub-item 2 tab content"
      ),
      tabItem("enrich",
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
              actionBttn("go", "GO"),
              tags$h3('Positive'),
              tags$h4("GO_Biological_Process_2021"),
              plotOutput("enrichment1"),
              tags$h4("GO_Molecular_Function_2021"),
              plotOutput("enrichment2"),
              tags$h4("Reactome_2022"),
              plotOutput("enrichment3"),
              tags$h4('Negative'),
              tags$h4("GO_Biological_Process_2021"),
              plotOutput("enrichment4"),
              tags$h4("KEGG_2021_Human"),
              plotOutput("enrichment5"),
              tags$h4("Reactome_2022"),
              plotOutput("enrichment6")
      )
    )
  )
  
  shinyApp(
    ui = dashboardPage(header, sidebar, body),
    server = function(input, output) {
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
      
      output$i<-renderPrint({
        input$i
      })
      
      gene_tissue<-reactive({
        origin_tissue = paste0(input$gene, input$origin_tissue)
      })
      
      col_scheme<-reactive({
        working_data<-working_data()
        col_scheme= rev(met.brewer('Austria', length(unique(working_data$tissue_2))))
        names(col_scheme) = unique(working_data$tissue_2)
      }) 
      
      output$topgene1<-renderPlot({
        working_data<-working_data()
        col_scheme= rev(met.brewer('Austria', length(unique(working_data$tissue_2))))
        names(col_scheme) = unique(working_data$tissue_2)
        max_gene_length = input$i
        max_gene_length<-as.numeric(max_gene_length)
        top_genes = working_data[1:max_gene_length,]
        top_genes$color = col_scheme[match(top_genes$tissue_2, names(col_scheme))]
        
        p<-ggplot(top_genes, aes(x=fct_reorder(gene_tissue_2, abs(bicor), .desc=T), y=bicor, fill=tissue_2)) + geom_col(position = position_dodge2(reverse = TRUE)) +  theme_classic() + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values = top_genes$color[order(abs(top_genes$bicor), decreasing = T)]) + xlab('') + ylab('bicor coefficent') + ggtitle('Top genes correlated with origin gene tissue')
        print(p)
      })
      
      observeEvent(
        input$bttn1, 
        {    # NAMED VECTOR CAUSES THE WARNING
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
              dplyr::mutate(freq = n / sum(n))%>%
              dplyr::arrange(desc(freq))
            
            
            tissue_freqs = binned_sig_prots %>% dplyr::group_by(qcat) %>% dplyr::summarise(sum(n))
            
            binned_sig_prots$tot_count = tissue_freqs$`sum(n)`[match(binned_sig_prots$qcat, tissue_freqs$qcat)]
            binned_sig_prots$qcat1 = paste0(binned_sig_prots$qcat, ' ', binned_sig_prots$tot_count, ' genes')
            col_scheme = rev(met.brewer('Austria', length(unique(sig_table$tissue_2))))
            names(col_scheme) = unique(sig_table$tissue)
            
            ggplot(binned_sig_prots, aes(x = "", y = freq, fill =tissue_2)) + 
              geom_bar(stat = "identity", width = 1, position = position_fill()) + blank_theme + theme(plot.title=element_text(size=25, face="bold")) +
              theme(axis.text.x=element_blank())+ scale_fill_manual(values=col_scheme) +
              coord_polar(theta = "y") + 
              facet_wrap( ~ qcat1)
            
            
          }) 
        })
      observeEvent(
        input$bttn1, 
        {    # NAMED VECTOR CAUSES THE WARNING
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
        })
      
      # pathway network
     
      
    
      
      # enrichment 
      enriched1 = reactive({
        annot <- working_data()

        select_tissue = input$tissue2
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
        
        
      })
      
      observeEvent(
        input$go, 
        {    # NAMED VECTOR CAUSES THE WARNING
          output$enrichment1<-renderPlot({
        enriched<-enriched1()
        plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
      })
      output$enrichment2<-renderPlot({
        enriched<-enriched1()
        plotEnrich(enriched[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
      })
      output$enrichment3<-renderPlot({
        enriched<-enriched1()
        plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
      })
        })
      
      
      enriched2 = reactive({
        annot <- working_data()
        
        select_tissue = input$tissue2
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
        
        
      })
      observeEvent(
        input$go, 
        {    # NAMED VECTOR CAUSES THE WARNING
          output$enrichment4<-renderPlot({
        enriched<-enriched2()
        plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
      })
      output$enrichment5<-renderPlot({
        enriched<-enriched2()
        plotEnrich(enriched[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
      })
      output$enrichment6<-renderPlot({
        enriched<-enriched2()
        plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
      })
        })
      

    }
  )
}

  
