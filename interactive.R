library(shiny)
library(dashboard)

ui <- dashboardPage(
  dashboardHeader(title = "My Dashboard"),
  dashboardSidebar(
    sidebarUserPanel("Let's start it"),
    sidebarMenu(
      # Setting id makes input$tabs give the tabName of currently-selected tab
      id = "tabs",
      menuItem("Instruction", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Settings",  tabName = "Settings",icon = icon("fa-solid fa-microscope")),
      menuItem("Let's see the pie chart", tabName = "Pie",icon = icon("fa-doutone fa-chart-pie")),
      menuItem("The pathway", tabName = "pathway",icon = icon("fa-doutone fa-diagram-project")),
      menuItem("The enrichemtn",tabName = "enrich", icon = icon("fa-doutone fa-chart-bar"))
    )
  ),
  dashboardBody(
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
                billboarderOutput("plot1")
              ),
              billboarderOutput("plot2"),
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
              tags$h4("GO_Biological_Process_2021"),
              plotOutput("enrichment1"),
              tags$h4("GO_Molecular_Function_2021"),
              plotOutput("enrichment2"),
              tags$h4("Reactome_2022"),
              plotOutput("enrichment3"),
              tags$h4("DSigDB"),
              plotOutput("enrichment4")
      ))))


server <- function(input, output) {
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
  
  
  output$plot1 = renderBillboarder({
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
    binned_sig_prots$tissue_2 <-  reorder(binned_sig_prots$tissue_2, binned_sig_prots$n, decreasing = T)
    #binned_sig_prots<-binned_sig_prots[!grepl(brain, binned_sig_prots$tissue_2)]
    billboarder(data = binned_sig_prots) %>% 
      bb_aes(tissue_2, freq) %>% 
      bb_piechart()%>% 
      bb_legend(position = 'right') %>%
      bb_color(palette = col_scheme)
    
    
  })
  output$plot2 = renderBillboarder({
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
    binned_sig_prots$tissue_2 <-  reorder(binned_sig_prots$tissue_2, binned_sig_prots$n, decreasing = T)
    

    
    billboarder(data = binned_sig_prots) %>% 
      bb_aes(tissue_2, freq) %>% 
      bb_piechart()%>% 
      bb_legend(position = 'right') %>%
      bb_color(palette = col_scheme)
  })
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


shinyApp(ui, server)
