library(shiny)
library(shinydashboard)
library(shinyalert)
library(dplyr)
library(purrr)
library(ggplot2)
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
library(qvalue)
library(RPostgres)
library(enrichR)
library(stats)
library(feather)
library(tidyr)
library(tibble)
library(gridExtra)
library(ggpubr)
library(writexl)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library("pathview")
library(cowplot)
library(europepmc)
library(fresh)
library(future)
library(progress)
allowWGCNAThreads()
## load data
load("mouse.RData")

## preset functions for analysis 
#generate pie chart with all tissues and genes
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
        e_chart(tissue_2, width='400px', height = NULL) %>%
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
#generate pie chart with origin tissue and gene removed
f2<-function(annot,origin_tissue, col_scheme){
  annot<-annot[!grepl(origin_tissue, annot$tissue_2),]
  sig_table=annot
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
        e_chart(tissue_2, width='400px', height = NULL) %>%
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

#  positive enrichment 
f_e1 <- function(annot,select_tissue,q){
  pp1 = annot[annot$qvalue<as.numeric(q),]
  pp1 = pp1[pp1$bicor>0,]
  pp1 = pp1[pp1$tissue_2 %in% select_tissue,]
  pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
  pp2 = pp1[1:pp1_length,]
  gg1 = pp2$gene_symbol_2
  
  setEnrichrSite("Enrichr")
  dbs <- listEnrichrDbs()
  dbs1 <- "GO_Biological_Process_2021"
  
  enriched <- enrichr(gg1, dbs1)
  
}

f_e2 <- function(annot,select_tissue,q){
  pp1 = annot[annot$qvalue<as.numeric(q),]
  pp1 = pp1[pp1$bicor>0,]
  pp1 = pp1[pp1$tissue_2 %in% select_tissue,]
  pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
  pp2 = pp1[1:pp1_length,]
  gg1 = pp2$gene_symbol_2
  
  setEnrichrSite("Enrichr")
  dbs <- listEnrichrDbs()
  dbs1 <- "Reactome_2022"
  
  enriched <- enrichr(gg1, dbs1)
  
  
}

# negative enrichment 
f_e3 <- function(annot,select_tissue,q){
  pp1 = annot[annot$qvalue<as.numeric(q),]
  pp1 = pp1[pp1$bicor<0,]
  pp1 = pp1[pp1$tissue_2 %in% select_tissue,]
  pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
  pp2 = pp1[1:pp1_length,]
  gg1 = pp2$gene_symbol_2
  
  setEnrichrSite("Enrichr")
  dbs <- listEnrichrDbs()
  dbs1 <- "GO_Biological_Process_2021"
  
  enriched <- enrichr(gg1, dbs1)
  
}

f_e4 <- function(annot,select_tissue,q){
  pp1 = annot[annot$qvalue<as.numeric(q),]
  pp1 = pp1[pp1$bicor<0,]
  pp1 = pp1[pp1$tissue_2 %in% select_tissue,]
  pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
  pp2 = pp1[1:pp1_length,]
  gg1 = pp2$gene_symbol_2
  
  setEnrichrSite("Enrichr")
  dbs <- listEnrichrDbs()
  dbs1 <- "Reactome_2022"
  
  enriched <- enrichr(gg1, dbs1)
}


# top genes
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