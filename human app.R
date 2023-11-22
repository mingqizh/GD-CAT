source('human global.R', local = TRUE)
options(future.globals.maxSize= 891289600)
## setting header & sidebar
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
    p("This is for human data."),
    radioButtons("Gender", "Choose Sex",c("Both"="both","Male"="male", "Female"="female")),
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
        Lung = "Lung", 
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
              width = 12,height="500px",
              tabPanel('All genes q<0.1', echarts4rOutput('plot.p1')),
              tabPanel('All genes q<0.01', echarts4rOutput('plot.p2')),
              tabPanel('All genes q<0.001', echarts4rOutput('plot.p3')),
              tabPanel('Origin-removed q<0.1', echarts4rOutput('plot.p4')),
              tabPanel('Origin-removed q<0.01', echarts4rOutput('plot.p5')),
              tabPanel('Origin-removed q<0.001', echarts4rOutput('plot.p6')),
              tabPanel('Table-1',DTOutput('table.t1'))
            ),
            shinydashboard::box(
              title = "Please click the pie chart body above to start the pathway section",
              width = 12,
              verbatimTextOutput('text'),
              #actionButton('btn', class = "btn-primary", 'Start Analysis'),
              #uiOutput('plots.en')
            ),
            shinydashboard::box(
              title = "Top GSEA-GO Activated and Suppressed",
              width = 12,
              sliderInput(inputId = "tp",label = "How many top pathways you want to see",value = 15, min = 1, max = 20),
              actionButton('dotb', class = "btn-primary", 'Start process'),
              plotOutput('dot',width  = "800px",height = "900px"),
              downloadButton("dotp", "Download Image"),
              plotOutput('nete',width  = "1200px",height = "900px"),
              downloadButton("ed", "Download Image")
            ),
            
            tabBox(title="Correlation with age or sex difference",width=12,height="600px",
                   tabPanel('Age and sex of the cohort',plotOutput('AS'),
                            downloadButton("ASP", "Download Image"))
            ),
            tabBox(title="Cell type",width=12,height="500px",
                   tabPanel('Cell type',echarts4rOutput('cell'))
            ),
            sliderInput(inputId = "topn",label = "How many top-ranked correlated genes do you want to see?",value = 30, min = 1, max = 50),
            tabBox(title="Top-N",width=12,height="500px",
                   tabPanel('Top-ranked genes',echarts4rOutput('plot.top1')),
                   tabPanel('Top-ranked genes without origin tissue',echarts4rOutput('plot.top2'))
            ),
            shinydashboard::box(
              title = "Scatter plot",
              width = 12,
              textInput("gene_sc","Input a gene for scatter plot, please use Official NBCI gene symbol",value = "ADIPOQ"),
              selectInput(
                "tissue_sc",
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
              actionButton('scb', class = "btn-primary", 'Start plot'),
              plotOutput('sc'),
              downloadButton("SCP", "Download Image")
            ),
            sliderInput(inputId = "within",label = "How many within tissue gene numbers you want to input in the network",value = 30, min = 1, max = 50),
            sliderInput(inputId = "external",label = "How many peripheral tissue gene numbers you want to input in the network",value = 100, min = 1, max = 150),
            shinydashboard::box(
              title = "Network Analysis",
              width = 12,
              actionButton('btn1', class = "btn-primary", 'Start Analysis'),
              tabPanel('Chart',plotOutput('net', width  = "1200px",height = "900px")),
              downloadButton("PDFPlot", "Download PDF"),
              downloadButton("netgene", "Download gene lists")
            )
          )
  )
}
# setting web body
body <- function(){
  dashboardBody(
    tabItems(
      t2()
    )
  )
}

#combine as UI
ui<-dashboardPage(header(), sidebar(), body())


# setting the server
server <- function(input, output, session) {
  observeEvent(input$tabs, {
    if(input$tabs=='t2'){
      shinyalert("Welcome to GD-CAT(Human)!", 
                 "If the page becomes unresponsive following a user action, it is likely a result of high website traffic. We'd appreciate your patience, as you are queued next in line.", 
                 type = "success")
    }
  })
  working_dataset<-eventReactive(input$import,{
    
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

    tryCatch({
      full_cors = bicorAndPvalue(tissue1, tissue2, use = 'p')
    }, error = function(e) {
      shinyalert("Oops!","Please check that you input the official NBCI gene symbol; another reason may be that no such gene is available in the dataset.", type = "error")
    })
      
    progress$set(value = 2)
    cor_table = reshape2::melt(full_cors$bicor)
    new_p = reshape2::melt(full_cors$p)
      
    colnames(cor_table) = c('gene_tissue_1', 'gene_tissue_2', 'bicor')
    #can drop here to clear CPU
    full_cors=NULL
    progress$set(value = 3)
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
      
    #sig_table åé¢éè¦ç¨å°
    sig_table = cor_table[cor_table$qvalue<0.1,]
    sig_table$qcat =ifelse(sig_table$qvalue<0.01, 'q<0.01', 'q<0.1')
    sig_table$qcat =ifelse(sig_table$qvalue<0.001, 'q<0.001', paste0(sig_table$qcat))
      
    progress$set(value = 4)
    
    
    
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
  
  AS<-reactive({
    working_dataset<-working_dataset()
    origin_gene<-input$origin_gene
    origin_tissue=input$origin_tissue
    origin_gene_tissue = paste0(origin_gene, '_', origin_tissue)
    
    met_dat = met_dat
    
    sex_table = met_dat
    sex_table$GTEx_ID = gsub('GTEX-', '', sex_table$SUBJID)
    sex_table$sexMF = ifelse(sex_table$SEX==1, 'M', 'F')
    new_trts = sex_table[sex_table$GTEx_ID %in% row.names(working_dataset),]
    table(new_trts$sexMF)
    #  F   M
    
    table(new_trts$AGE)
    age_factors = as.data.frame(table(new_trts$AGE))
    age_factors$seq_num = seq(length(row.names(age_factors)))
    new_trts$age_num = age_factors$seq_num[match(new_trts$AGE, age_factors$Var1)]
    
    trait_df = as.data.frame(working_dataset[,colnames(working_dataset) %in% origin_gene_tissue])
    
    colnames(trait_df) = 'gene_tissue'
    row.names(trait_df) = row.names(working_dataset)
    trait_df$age = new_trts$AGE[match(row.names(trait_df), new_trts$GTEx_ID)]
    trait_df$age_n = new_trts$age_num[match(row.names(trait_df), new_trts$GTEx_ID)]
    
    trait_df$sex = new_trts$sexMF[match(row.names(trait_df), new_trts$GTEx_ID)]
    new_cors1 = bicorAndPvalue(trait_df$gene_tissue, trait_df$age_n, use = 'p')
    age_plot1 = ggplot(trait_df, aes(x=factor(age_n), y=gene_tissue)) + geom_point() + theme_classic() + geom_smooth(aes(x=age_n, y=gene_tissue),method = 'lm')+ xlab('Age group (years)') + ylab(paste0('Normalized ', origin_gene_tissue, ' expression')) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
      scale_x_discrete(labels = c('20-29', '30-39', '40-49', '50-59', '60-69', '70-79'))  +  ggtitle(paste0(origin_gene_tissue, ' ~ Age, bicor = ', round(new_cors1$bicor, 3), 'p = ', signif(new_cors1$p, 3)))
    sex_plot1 = ggplot(trait_df, aes(x=sex, y=gene_tissue, fill=sex)) +
      geom_violin(trim=FALSE )+
      geom_boxplot(width=0.1) + theme_minimal() + ylab(paste0('Normalized ', origin_gene_tissue, ' expression')) + xlab('Reported sex') + theme(legend.position = "none")+ scale_fill_manual(values=c('darkslategray1',"gold")) + stat_compare_means()
    final_plot_metadata = grid.arrange(age_plot1, sex_plot1, ncol=2)
    
  })
  output$AS<-renderPlot({
    AS()
  })
  
  output$ASP <- downloadHandler(
    filename = function() {
      "plot.png"
    },
    content = function(file) {
      ggsave(file, plot = AS(), device = "png", width = 9, height = 6)
    }
  )
  
  
  #cell type
  output$cell<- renderEcharts4r({
    top_genes1<-get_cell(working_dataset(),sig_table(), input$origin_gene, input$origin_tissue, col_scheme())
    
    top_genes1 %>%
      arrange(desc(abs(bicor))) %>%
      e_chart(cell_tissue, width='300px', height = NULL) %>%
      e_bar(bicor) %>%
      e_legend(show=F) %>%
      e_add_nested("itemStyle", color) %>%
      e_grid(bottom="150px") %>%
      e_x_axis(axisLabel = list(interval = 0, rotate = 45)) %>%
      e_y_axis(name = "Bicor Values")%>%
      e_datazoom(type='inside') %>%
      e_tooltip(
        trigger = 'item',
        axisPointer = list(
          type = "shadow",
          axis='x'
        )
      )%>%
      e_toolbox()%>%
      e_toolbox_feature(feature = "saveAsImage",title='Save')
  }
  )
  
  # top n
  output$plot.top1<- renderEcharts4r({
    top_genes1<-get_top_genes1(sig_table(),input$topn, col_scheme())
    
    top_genes1 %>%
      arrange(desc(abs(bicor))) %>%
      e_chart(gene_tissue_2, width='300px', height = NULL) %>%
      e_bar(bicor) %>%
      e_legend(show=F) %>%
      e_add_nested("itemStyle", color) %>%
      e_grid(bottom="150px") %>%
      e_x_axis(axisLabel = list(interval = 0, rotate = 45)) %>%
      e_y_axis(name = "Bicor Values")%>%
      e_datazoom(type='inside') %>%
      e_tooltip(
        trigger = 'item',
        axisPointer = list(
          type = "shadow",
          axis='x'
        )
      )%>%
      e_toolbox()%>%
      e_toolbox_feature(feature = "saveAsImage",title='Save')%>%
      e_on(
        list(seriesName = bicor),
        "function(x){
          //alert(Object.keys(x));
          var msg = [x.seriesIndex,x.seriesName,x.name,x.dataIndex]
          //alert(msg)
          Shiny.setInputValue('selected_g',x.name, {priority: 'event'})
          Shiny.setInputValue('selected_gene',x.seriesName, {priority: 'event'});
        }"
      )
  }
  )
  output$gene1 <- renderText({
    print(c('You selected gene-tissue:',input$selected_gene))
  })
  output$plot.top2<- renderEcharts4r({
    top_genes2<-get_top_genes2(sig_table(),input$topn,input$origin_tissue, col_scheme())
    
    top_genes2 %>%
      arrange(desc(abs(bicor))) %>%
      e_chart(gene_tissue_2, width='300px', height = NULL) %>%
      e_bar(bicor) %>%
      e_legend(show=F) %>%
      e_add_nested("itemStyle", color) %>%
      e_grid(bottom="100px") %>%
      e_x_axis(axisLabel = list(interval = 0, rotate = 45)) %>%
      e_y_axis(name = "Bicor Values")%>%
      e_datazoom(type='inside') %>%
      e_tooltip(
        trigger = 'item',
        axisPointer = list(
          type = "shadow",
          axis='x'
        )
      )%>%
      e_toolbox()%>%
      e_toolbox_feature(feature = "saveAsImage",title='Save')
  }
  )
  
  # scatter plot
  observeEvent( input$scb, {
    if(!is.null(input$gene_sc)){
      
      isolate({
        origin_gene_tissue = paste0(input$origin_gene, '_', input$origin_tissue)
        select_gene2 = paste0(input$gene_sc, '_', input$tissue_sc)
        
      })
      
      working_dataset<-working_dataset()
      scat_gene_set1 = working_dataset[,  colnames(working_dataset) %in% origin_gene_tissue]
      scat_gene_set2 = working_dataset[,  colnames(working_dataset) %in% select_gene2]
      
      scat_df = as.data.frame(cbind(scat_gene_set1, scat_gene_set2))
      tryCatch({
        colnames(scat_df) = c('sgene_1', 'sgene2')
        cc1 = bicorAndPvalue(scat_df$sgene2, scat_df$sgene_1)
        sc<-ggplot(scat_df, aes(x=sgene_1, y=sgene2)) + theme_classic() +  geom_hex(bins = 100) + scale_fill_distiller(palette = "Blues", direction=-1)+ geom_smooth(method = 'lm', color = "darkorange")+  ggtitle(paste0(origin_gene_tissue, ' ~ ', select_gene2, ' r=', round(cc1$bicor, 3), ' p=', signif(cc1$p, 3))) + xlab(paste0('Normalized ', origin_gene_tissue, ' expression')) + ylab(paste0('Normalized ', select_gene2, ' expression')) + theme(aspect.ratio=1)
        output$sc<-renderPlot({
          sc
        }
        )
        output$SCP <- downloadHandler(
          filename = function() {
            "Scatter.png"
          },
          content = function(file) {
            ggsave(file, plot = sc, device = "png", width = 9, height = 6)
          }
        )
        
      }, error= function(e) {
        shinyalert("Oops!","Plese input a official NBCI gene symbol.", type = "error")
      })
    } else {
      
      shinyalert("Warning!", "Please first input an gene-tissue.", type = "warning")
    }
  })
  
  observeEvent(input$btn1,{
    progress <- Progress$new(session, min=0, max=5)
    on.exit(progress$close())
    progress$set(message = 'Gnerating the network plot',
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
    progress$set(value = 2)
    
    origin_pull = sig_table[sig_table$tissue_2==origin_tissue,]
    orig_network_genes = as.vector(origin_pull$gene_tissue_2[1:number_orig_gene])
    peripheral_pull = sig_table[!sig_table$tissue_2==origin_tissue,]
    periph_network_genes = as.vector(peripheral_pull$gene_tissue_2[1:number_peripheral_genes])
    sql_pull_list = c(origin_gene_tissue, orig_network_genes, periph_network_genes)
      
      
      
    tissue1 <- working_dataset[, colnames(working_dataset)  %in% sql_pull_list]
    full_cors = bicorAndPvalue(tissue1, tissue1, use = 'p')
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
    cor_table$gene_symbol_2 = gsub("\\_.*","",cor_table$gene_tissue_2)
    cor_table$tissue_2 = gsub(".*_","",cor_table$gene_tissue_2)
    cor_table = cor_table[!is.na(cor_table$tissue_2),]
    sql_pull_2 = cor_table
      
    
    progress$set(value = 3)
    ######################################
    #the new SQL pull table is not listed as sql_pull_2
    #this will take work on your end to interface with pulling the correct genes from tissue_1 based on the number selected in the app.  I arbitrarily chose 10 here and you obviously wouldn't jsut apply top_n since it will differ.  The final list is called network_gene_lsit
    map1 = sql_pull_2
    map1$cols = col_scheme[match(map1$tissue_1, names(col_scheme))]
    map2 = dcast(map1, gene_tissue_1 ~ gene_tissue_2, value.var = 'bicor', fun.aggregate = mean)
    row.names(map2) = map2$gene_tissue_1
    map3<-map2
    
    map2$gene_tissue_1 = NULL
    cols_set = as.data.frame(row.names(map2))
    row.names(cols_set) = row.names(map2)
    cols_set$cols = map1$cols[match(row.names(cols_set), map1$gene_tissue_1)]
    cols_set$cols = ifelse(row.names(cols_set) %in% origin_gene_tissue, 'gray5', paste0(cols_set$cols))
    progress$set(value = 4)
    output$net<-renderPlot({
      qgraph(map2, minimum = 0.5, cut = 0.85, vsize = 3, color=cols_set$cols, legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=4, directed=F, labels = colnames(map2)) 
    })
    output$PDFPlot <- downloadHandler(
      filename = function() {
        paste("Net-",origin_gene_tissue,"-", input$Gender,"-",  Sys.Date(), ".pdf", sep="")
      },
      content = function(file) {
        pdf(file, width = 20, height = 20)
        plot(qgraph(map2, minimum = 0.5, cut = 0.85, vsize = 3, color=cols_set$cols, legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=4, directed=F, labels = colnames(map2)) )
        dev.off()
      }
    )
    output$netgene<- downloadHandler(
      filename = function() {
        paste("Gene list for network plot-",origin_gene_tissue,"-", input$Gender,"-", Sys.Date(), ".xlsx", sep="")
      },
      content = function(file) {
        write_xlsx(map3, path = file)
      }
    )
  })
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
    pie2()[[3]]
  })
  output$plot.p5 <- renderEcharts4r({
    pie2()[[2]]
  })
  output$plot.p6 <- renderEcharts4r({
    pie2()[[1]]
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
                   lengthMenu = list(c(10,25,50),
                                     c(10,25,50,"All")) ))
  
  observeEvent(input$dotb,{
    if(!is.null(input$selected_tissue)){
      progress <- Progress$new(session, min=0, max=5)
      on.exit(progress$close())
      progress$set(message = 'Processing the object for pathways',
                   detail = 'It will take around 2-5 minutes depending on usage')
      progress$set(value = 1)
      isolate({
        sig_table<-sig_table()
        select_tissue = input$selected_tissue
        pie_bin = as.numeric(input$selected_q)
      })
      binned_sig_prots= sig_table %>%
        dplyr::group_by(qcat, tissue_2) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::mutate(freq = n / sum(n))
      
      tissue_freqs = binned_sig_prots %>% dplyr::group_by(qcat) %>% dplyr::summarise(sum(n))
      
      binned_sig_prots$tot_count = tissue_freqs$`sum(n)`[match(binned_sig_prots$qcat, tissue_freqs$qcat)]
      binned_sig_prots$qcat1 = paste0(binned_sig_prots$qcat, ' ', binned_sig_prots$tot_count, ' genes')
      progress$set(value = 2)
      tissue_table = binned_sig_prots[binned_sig_prots$qcat =='q<0.1',]
      
      
      pp1 = sig_table[sig_table$tissue_2 %in% select_tissue,]
      pp1 = pp1[pp1$pvalue < pie_bin,]
      #fc_dko = scales::rescale(pp1$bicor, to=c(-3, 3))
      fc_dko = pp1$bicor
      progress$set(value = 3)
      ## match each fold change value with the corresponding gene symbol
      names(fc_dko) <- pp1$gene_symbol_2
      #Next we need to order the fold changes in decreasing order. To do this we'll use the sort() function, which takes a vector as input. This is in contrast to Tidyverse's arrange(), which requires a data frame.
        
      ## Sort fold changes in decreasing order
      fc_dko <- sort(fc_dko, decreasing = TRUE)
        
      organism = "org.Hs.eg.db"
        
        ## Org.Hs.eg.db https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html
        
      gse <-gseGO(
        geneList=fc_dko,
        ont = "ALL",
        OrgDb= organism,
        keyType = "SYMBOL",
        exponent = 1,
        minGSSize = 2,
        maxGSSize = 500,
        eps = 0,
        pvalueCutoff = 1,
        pAdjustMethod = "BH") 
      str(gse)
        
      progress$set(value = 4)
      isolate({
        origin_gene<-input$origin_gene
        origin_tissue=input$origin_tissue
        origin_gene_tissue = paste0(origin_gene, '_', origin_tissue)
        number<-as.numeric(input$tp)
      })
      
      output$dot<-renderPlot({
        dotplot(gse, showCategory=number, split=".sign", color = "pvalue") + facet_grid(.~.sign) +
          ggtitle(paste0('GSEA pathways from positive and negative correlations ',  origin_gene_tissue, ' in ', input$selected_tissue, ' ', input$Gender))
      })
      
      output$dotp <- downloadHandler(
        filename = function() {
          paste("GSEA Pathway-",origin_gene_tissue,' in ', select_tissue, ' ', input$Gender,"-", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, width = 10, height = 10)
          plot(dotplot(gse, showCategory=number, split=".sign", color = "pvalue") + facet_grid(.~.sign) +
                 ggtitle(paste0('GSEA pathways from positive and negative correlations ',  origin_gene_tissue, ' in ', select_tissue, ' ', input$Gender))
          )
          dev.off()
        }
      )
      x2<- pairwise_termsim(gse)
      output$nete<-renderPlot({
        emapplot(x2, showCategory=number, color = "pvalue")+ ggtitle("Relationship between the top most significantly GSE - GO terms (padj.)")
      })
      progress$set(value = 5)
      output$ed <- downloadHandler(
        filename = function() {
          paste("GSEA Network-",origin_gene_tissue,' in ', select_tissue, ' ', input$Gender,"-", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, width = 10, height = 10)
          plot(emapplot(x2, showCategory=number, color = "pvalue")+ ggtitle("Relationship between the top most significantly GSE - GO terms (padj.)")
          )
          dev.off()
        }
      )
    }else {
      shinyalert("Warning!", "Please first click the Pie chart body to selected an tissue.", type = "warning")
    }
  })
  
  # click
  output$text <- renderText({
    print(c('You selected tissue:',input$selected_tissue,', and q<',input$selected_q))
  })
  # enrichment
  plots<-reactiveVal()
  en1<-reactiveVal()
  en2<-reactiveVal()
  en3<-reactiveVal()
  en4<-reactiveVal()
  observeEvent( input$btn, {
    if(!is.null(input$selected_tissue)){
      progress <- Progress$new(session, min=0, max=5)
      on.exit(progress$close())
      progress$set(message = 'Generating the enrichements',
                   detail = 'Just wait a second')
      progress$set(value = 1)
      tryCatch({
        # Your plot code here
        isolate({
          enriched1<-f_e1(sig_table(),input$selected_tissue,input$selected_q)
        })
        
        plots_en1<- enriched1 %>%
          purrr::map(~plotEnrich(.x, showTerms = 10, numChar = 30, y = "Count", orderBy = "P.value") 
                     + ggtitle(paste0('Positive gene correlations with ', input$origin_gene, ' ',input$origin_tissue, ' ', 'GO_Biological_Process_2021')))
        progress$set(value = 2)
        en1(enriched1)
        isolate({
          enriched2<-f_e2(sig_table(),input$selected_tissue,input$selected_q)
        })
        
        plots_en2<- enriched2 %>%
          purrr::map(~plotEnrich(.x, showTerms = 10, numChar = 30, y = "Count", orderBy = "P.value")
                     + ggtitle(paste0('Positive gene correlations with ', input$origin_gene, ' ', input$origin_tissue, ' ','Reactome_2022')))
        en2(enriched2)
        isolate({
          enriched3<-f_e3(sig_table(),input$selected_tissue,input$selected_q)
        })
        
        plots_en3<- enriched3 %>%
          purrr::map(~plotEnrich(.x, showTerms = 10, numChar = 30, y = "Count", orderBy = "P.value") 
                     + ggtitle(paste0('Negative gene correlations with ', input$origin_gene, ' ',input$origin_tissue, ' ', 'GO_Biological_Process_2021')))
        
        en3(enriched3)
        isolate({
          enriched4<-f_e4(sig_table(),input$selected_tissue,input$selected_q)
        })
        
        plots_en4<- enriched4 %>%
          purrr::map(~plotEnrich(.x, showTerms = 10, numChar = 30, y = "Count", orderBy = "P.value")
                     + ggtitle(paste0('Negative gene correlations with ', input$origin_gene, ' ', input$origin_tissue, ' ','Reactome_2022')))
        progress$set(value = 3)
        en4(enriched4)
        plots_all<-c(plots_en1,plots_en2, plots_en3,plots_en4)
        plots(plots_all)
        output[['plots.en']]<-renderUI({
          plot_output_list <- lapply(1:length(plots_all), function(i) {
            plotname <- paste("en", i, sep="")
            plotOutput(plotname)
          })
          plot_output_list$btn_down<- downloadButton("download", "Download Image")
          plot_output_list$table<-downloadButton("table", "Download Table")
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
        
        shinyalert("Oops","The selected tissue do not contain enough genes to generate the negative enrichment. Please select again.", type = "error")
        
        tryCatch({
          enriched1<-f_e1(sig_table(),input$selected_tissue,input$selected_q)
          plots_en1<- enriched1 %>%
            purrr::map(~plotEnrich(.x, showTerms = 10, numChar = 30, y = "Count", orderBy = "P.value") 
                       + ggtitle(paste0('Positive gene correlations with ', input$origin_gene, ' ',input$origin_tissue, ' ', 'GO_Biological_Process_2021')))
          progress$set(value = 2)
          en1(enriched1)
          enriched2<-f_e2(sig_table(),input$selected_tissue,input$selected_q)
          plots_en2<- enriched2 %>%
            purrr::map(~plotEnrich(.x, showTerms = 10, numChar = 30, y = "Count", orderBy = "P.value")
                       + ggtitle(paste0('Positive gene correlations with ', input$origin_gene, ' ', input$origin_tissue, ' ','Reactome_2022')))
          en2(enriched2)
          plots_all<-c(plots_en1,plots_en2)
          plots(plots_all)
          output[['plots.en']]<-renderUI({
            plot_output_list <- lapply(1:length(plots_all), function(i) {
              plotname <- paste("en", i, sep="")
              plotOutput(plotname)
            })
            plot_output_list$btn_down<- downloadButton("download", "Download Image")
            plot_output_list$table<-downloadButton("table", "Download Table")
            do.call(tagList, plot_output_list)
          }
          )
          
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
          shinyalert("Oops!","The selected tissue do not contain enough genes to generate positive the enrichment. Please select again.", type = "error")
        })
        
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
        ggsave( paste0('plot',i,'.png'), plot = plots[[i]], device = "png", width = 12, height = 9)
      }
      zip::zip(file,paste0('plot',1:length(plots),'.png'))
    }
  )
  
  output$table<- downloadHandler(
    filename = function() {
      paste("Enrichments-", Sys.Date(), ".xlsx", sep="")
    },
    content = function(file) {
      write_xlsx(list(positive_GO_Biological_process = en1()[[1]], Positive_Reactome = en2()[[1]], Negative_GO_Biological_process = en3()[[1]], Negative_Reactome = en4()[[1]]), path = file)
    }
  )
  #tip
  observeEvent(input$tabs, {
    if(input$tabs=='t2'){
      showNotification("Click the legend (on the right of the chart body) of the Pie chart to toggle the display of the series.",duration = NULL,type="message")
    }
  })
  
}

#run the app
shinyApp(ui, server)

