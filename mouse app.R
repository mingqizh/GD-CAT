source("mouse global.R", local = TRUE)

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
    p("This is for mouse data."),
    radioButtons("diet", "Choose Diet",c("HF"="HF","Chow"="Chow")),
    textInput("origin_gene","origin_gene, Official NBCI gene symbol",value = "Adipoq"),
    conditionalPanel(
      condition = "input.diet=='HF'",
      selectInput(
        "origin_tissue",
        "Select a tissue",
        c(Adipose="adipose",
          Aorta = "aorta", 
          Hypothlamus = "hypothlamus", 
          Heart = "heart",
          Liver = "liver",
          Intestine="intestine",
          Muscle="muscle")
      )
    ),
    conditionalPanel(
      condition = "input.diet=='Chow'",
      selectInput(
        "origin_tissue",
        "Select a tissue",
        c(Adipose="adipose",
          Aorta = "aorta", 
          Bone = "bone", 
          Heart = "heart",
          Liver = "liver")
      )
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
              title = "Please click the pie chart body above to start the pathway section",
              width = 12,
              verbatimTextOutput('text')
              #actionButton('btn', class = "btn-primary", 'Start Analysis'),
              #uiOutput('plots.en')
            ),
            shinydashboard::box(
              title = "Top GSEA-GO Activated and Suppressed",
              width = 12,
              sliderInput(inputId = "tp",label = "How many top pathways you want to see",value = 15, min = 1, max = 20),
              actionButton('dotb', class = "btn-primary", 'Start process'),
              plotOutput('dot',width  = "900px",height = "900px"),
              downloadButton("dotp", "Download Image"),
              plotOutput('nete',width  = "900px",height = "900px"),
              downloadButton("ed", "Download Image")
            ),
            sliderInput(inputId = "topn",label = "How many top-ranked correlated genes do you want to see?",value = 30, min = 1, max = 50),
            tabBox(title="Top-N",width=12,height="500px",
                   tabPanel('Top-ranked genes',echarts4rOutput('plot.top1')),
                   tabPanel('Top-ranked genes without origin tissue',echarts4rOutput('plot.top2'))
            ),
            
            
            tabBox(
              title = "Clinic trait",
              width = 12,
              tabPanel('Trait Table',DTOutput('trait_t'))
            ),
            sliderInput(inputId = "trait",label = "How many traits you want to see for the correlations",value = 30, min = 1, max = 50),
            shinydashboard::box(
              title = "Trait correlations",
              width = 12,
              actionButton('Mbtn1', class = "btn-primary", 'Start Analysis'),
              tabPanel('Correlations',plotOutput('trait_c', width  = "800px",height = "1000px")),
              downloadButton("trait_b", "Download PDF")
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
ui<-dashboardPage(skin="green",header(), sidebar(), body())


# setting the server
server <- function(input, output, session) {
  observeEvent(input$tabs, {
    if(input$tabs=='t2'){
      shinyalert("Welcome to GD-CAT(Mouse)!", 
                 "If the page becomes unresponsive following a user action, it is likely a result of high website traffic. We kindly request your patience, as you are queued next in line.", 
                 type = "success")
    }
  })
  working_dataset<-eventReactive(input$import,{
    
    isolate({
      if(input$diet == "HF"){
        working_dataset<-HF
      }else{
        working_dataset<-Chow
      }
      working_dataset
    })
  })
  
  working_traits<-eventReactive(input$import,{
    
    isolate({
      if(input$diet == "HF"){
        working_traits<-HT
      }else{
        working_traits<-CT
      }
      working_traits
    })
  })
  plan(multisession)
  
  sig_table <- eventReactive(input$import, {
    # Create a promise to encapsulate the long-running computation
    progress <- Progress$new(session, min=0, max=5)
    on.exit(progress$close())
    progress$set(message = 'Pre-processing raw data',
                 detail = 'It will take around 2-5 minutes depending on usage')
    progress$set(value = 1)
    isolate({
      origin_gene <- input$origin_gene
      origin_tissue <- input$origin_tissue
      origin_gene_tissue <- paste0(origin_gene, '_', origin_tissue)
      working_dataset <- working_dataset()
    })
    promise <- future({
      tissue2 <- working_dataset
      tissue1 <- working_dataset[, colnames(working_dataset) %in% origin_gene_tissue]
      
      tryCatch({
        full_cors <- bicorAndPvalue(tissue1, tissue2, use = 'p')
      }, error = function(e) {
        shinyalert("Oops!", "Please check that you input the official NCBI gene symbol; another reason may be that no such gene is available in the dataset.", type = "error")
      })
      
      cor_table <- reshape2::melt(full_cors$bicor)
      new_p <- reshape2::melt(full_cors$p)
      colnames(cor_table) <- c('gene_tissue_1', 'gene_tissue_2', 'bicor')
      full_cors <- NULL
      
      cor_table$pvalue <- signif(new_p$value, 3)
      cor_table$bicor <- round(cor_table$bicor, 3)
      cor_table <- cor_table[cor_table$gene_tissue_1 != cor_table$gene_tissue_2, ]
      cor_table$qvalue <- signif(p.adjust(cor_table$pvalue, "BH"), 3)
      cor_table <- cor_table[order(cor_table$qvalue, decreasing = FALSE), ]
      cor_table <- na.omit(cor_table)
      cor_table$gene_symbol_1 <- gsub("\\_.*", "", cor_table$gene_tissue_1)
      cor_table$tissue_1 <- gsub(".*_", "", cor_table$gene_tissue_1)
      cor_table <- cor_table[!is.na(cor_table$tissue_1), ]
      
      cor_table$gene_symbol_2 <- gsub("\\_.*", "", cor_table$gene_tissue_2)
      cor_table$tissue_2 <- gsub(".*_", "", cor_table$gene_tissue_2)
      cor_table <- cor_table[!is.na(cor_table$tissue_2), ]
      # sig_table needed for later use
      sig_table <- cor_table[cor_table$qvalue < 0.1, ]
      sig_table$qcat <- ifelse(sig_table$qvalue < 0.01, 'q<0.01', 'q<0.1')
      sig_table$qcat <- ifelse(sig_table$qvalue < 0.001, 'q<0.001', paste0(sig_table$qcat))
      
      return(sig_table)
    })
    progress$set(value = 2)
    # Execute the promise and retrieve the result using future's value()
    sig_table <- future::value(promise)
    progress$set(value = 3)
    # Return the result
    return(sig_table)
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
  
  
  # top n
  observeEvent(input$import,{
    isolate({
      origin_gene = input$origin_gene
      origin_tissue = input$origin_tissue
      origin_gene_tissue = paste0(origin_gene, '_', origin_tissue)
      working_dataset<-working_dataset()
    })
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
  })
  
  
  # clinical table
  
  cor_table1<-reactive({
    working_traits<-working_traits()
    origin_gene = input$origin_gene
    origin_tissue = input$origin_tissue
    origin_gene_tissue = paste0(origin_gene, '_', origin_tissue)
    working_dataset<-working_dataset()
    tissue1 <- working_dataset[,colnames(working_dataset) %in% origin_gene_tissue]
    full_cors1 = bicorAndPvalue(tissue1, working_traits, use = 'p')
    cor_table1 = reshape2::melt(full_cors1$bicor)
    new_p = reshape2::melt(full_cors1$p)
    cor_table1$Var1=NULL
    colnames(cor_table1) = c( 'trait_name', 'bicor')
    #can drop here to clear CPU
    full_cors1=NULL
    
    cor_table1$pvalue = signif(new_p$value, 3)
    cor_table1 = cor_table1[order(cor_table1$pvalue, decreasing = F),]
    cor_table1
  })
  
  output$trait_t <-renderDT(
    cor_table1(),
    extensions = 'Buttons',
    options = list(dom = 'Blfrtip',
                   scrollX=TRUE,
                   autoWidth=FALSE,
                   buttons = c('copy','csv','excel'),
                   lengthMenu = list(c(10,25,50),
                                     c(10,25,50,"All")) ))
  observeEvent(input$Mbtn1,{
    isolate({
      topn<-input$trait
      origin_gene = input$origin_gene
      origin_tissue = input$origin_tissue
      origin_gene_tissue = paste0(origin_gene, '_', origin_tissue)
      working_dataset<-working_dataset()
    })
    cor_table1<-cor_table1()
    cc1 = cor_table1[1:topn,]
    cc1$direction = ifelse(cc1$bicor>0, 'positive correlation', 'negative correlation')
    cc1$logp = -log10(cc1$pvalue)
    output$trait_c<-renderPlot({
      ggdotchart(cc1, x = "trait_name", y = "logp",
                 color = "direction",                                # Color by groups
                 palette = c(  "#FC4E07", "#00AFBB"), # Custom color palette
                 sorting = "descending",                       # Sort value in descending order
                 add = "segments",                             # Add segments from y = 0 to dots
                 rotate = TRUE,                                # Rotate vertically
                 group = "direction",                                # Order by groups
                 dot.size = 6,                                 # Large dot size
                 label = round(cc1$bicor, 2),                        # Add mpg values as dot labels
                 font.label = list(color = "black", size = 9, 
                                   vjust = 0.5),               # Adjust label parameters
                 ggtheme = theme_pubr()                        # ggplot2 theme
      ) + ylab('-log10(pvalue)') + xlab('') + ggtitle(paste0(origin_gene_tissue, ' ~ trait correlations',' ', input$diet))
    })
    
    output$trait_b <- downloadHandler(
      filename = function() {
        paste("Correlation-",origin_gene_tissue,"-", input$diet," diet-", Sys.Date(), ".pdf", sep="")
      },
      content = function(file) {
        pdf(file)
        plot(ggdotchart(cc1, x = "trait_name", y = "logp",
                        color = "direction",                                # Color by groups
                        palette = c(  "#FC4E07", "#00AFBB"), # Custom color palette
                        sorting = "descending",                       # Sort value in descending order
                        add = "segments",                             # Add segments from y = 0 to dots
                        rotate = TRUE,                                # Rotate vertically
                        group = "direction",                                # Order by groups
                        dot.size = 6,                                 # Large dot size
                        label = round(cc1$bicor, 2),                        # Add mpg values as dot labels
                        font.label = list(color = "black", size = 9, 
                                          vjust = 0.5),               # Adjust label parameters
                        ggtheme = theme_pubr()                        # ggplot2 theme
        ) + ylab('-log10(pvalue)') + xlab('') + ggtitle(paste0(origin_gene_tissue, ' ~ trait correlations',' ', input$diet)))
        dev.off()
      }
    )
    
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
    progress$set(value = 2)
    
    
    origin_gene_tissue = paste0(origin_gene, '_', origin_tissue)
    col_scheme<- col_scheme()
    promise<-future({
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
      return(sql_pull_2)
    })
    sql_pull_2<-future::value(promise)
    
    
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
    
    progress$set(value = 3)
    output$net<-renderPlot({
      qgraph(map2, minimum = 0.5, cut = 0.85, vsize = 3, color=cols_set$cols, legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=4, directed=F, labels = colnames(map2)) 
    })
    output$PDFPlot <- downloadHandler(
      filename = function() {
        paste("Net-", origin_gene_tissue,"-", input$diet," diet-",  Sys.Date(), ".pdf", sep="")
      },
      content = function(file) {
        pdf(file, width = 20, height = 20)
        plot(qgraph(map2, minimum = 0.5, cut = 0.85, vsize = 3, color=cols_set$cols, legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=4, directed=F, labels = colnames(map2)) )
        dev.off()
      }
    )
    output$netgene<- downloadHandler(
      filename = function() {
        paste("Gene list for network plot-" ,origin_gene_tissue,"-", input$diet," diet-", Sys.Date(), ".xlsx", sep="")
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
    options = list(
      dom = 'Blfrtip',
      scrollX = TRUE,
      autoWidth = FALSE,
      buttons = list(
        list(
          extend = 'copy',
          exportOptions = list(modifier = list(page = 'all'))
        ),
        list(
          extend = 'csv',
          exportOptions = list(modifier = list(page = 'all'))
        ),
        list(
          extend = 'excel',
          exportOptions = list(modifier = list(page = 'all'))
        )
      ),
      lengthMenu = list(c(10, 25, 50, -1), c(10, 25, 50, "All"))
    ),
    server = FALSE
  )
  observeEvent(input$dotb,{
    if(!is.null(input$selected_tissue)){
      progress <- Progress$new(session, min=0, max=5)
      on.exit(progress$close())
      progress$set(message = 'Processing the object for GSEA pathways',
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
      
      tissue_table = binned_sig_prots[binned_sig_prots$qcat =='q<0.1',]
      
      
      pp1 = sig_table[sig_table$tissue_2 %in% select_tissue,]
      pp1 = pp1[pp1$pvalue < pie_bin,]
      head(pp1)
      #fc_dko = scales::rescale(pp1$bicor, to=c(-3, 3))
      fc_dko = pp1$bicor
      
      ## match each fold change value with the corresponding gene symbol
      names(fc_dko) <- pp1$gene_symbol_2
      progress$set(value = 2)
      #Next we need to order the fold changes in decreasing order. To do this we'll use the sort() function, which takes a vector as input. This is in contrast to Tidyverse's arrange(), which requires a data frame.
      promise<-future({
        ## Sort fold changes in decreasing order
        fc_dko <- sort(fc_dko, decreasing = TRUE)
        
        organism = "org.Mm.eg.db"
        
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
        return(gse)
        
      })
      gse<-future::value(promise)
      
      isolate({
        origin_gene<-input$origin_gene
        origin_tissue=input$origin_tissue
        origin_gene_tissue = paste0(origin_gene, '_', origin_tissue)
        number<-as.numeric(input$tp)
      })
      output$dot<-renderPlot({
        
        
        dotplot(gse, showCategory=number, split=".sign", color = "pvalue") + facet_grid(.~.sign) +
          ggtitle(paste0('GSEA pathways from positive and negative correlations ',  origin_gene_tissue, ' in ', input$selected_tissue,' ',  input$diet))
      })
      progress$set(value = 3)
      output$dotp <- downloadHandler(
        filename = function() {
          paste("GSEA Pathway-",origin_gene_tissue,' in ', select_tissue, ' ', input$diet," diet-",  Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, width = 10, height = 10)
          plot(dotplot(gse, showCategory=number, split=".sign", color = "pvalue") + facet_grid(.~.sign) +
                 ggtitle(paste0('GSEA pathways from positive and negative correlations ',  origin_gene_tissue, ' in ', select_tissue, ' ', input$diet))
          )
          dev.off()
        }
      )
      x2<- pairwise_termsim(gse)
      output$nete<-renderPlot({
        emapplot(x2, showCategory = number, color = "pvalue")+ ggtitle("Relationship between the top most significantly GSE - GO terms (padj.)")
      })
      progress$set(value = 4)
      output$ed <- downloadHandler(
        filename = function() {
          paste("GSEA Network-", origin_gene_tissue,' in ', select_tissue, ' ', input$diet," diet-", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file, width = 10, height = 10)
          plot(emapplot(x2, showCategory = number, color = "pvalue")+ ggtitle("Relationship between the top most significantly GSE - GO terms (padj.)")
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
  #tip
  observeEvent(input$tabs, {
    if(input$tabs=='t2'){
      showNotification("Click the legend (on the right of the chart body) of the Pie chart to toggle the display of the series.",duration = NULL,type="message")
    }
  })
  
}

#run the app
shinyApp(ui, server)

