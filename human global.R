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
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library("pathview")
library(cowplot)
library(future)
library(progress)
allowWGCNAThreads()
## load data
load("working_dataset.RData")
met_dat = read.delim('GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
sc_matrix = read.csv('full pan-tissue DECON abundances.csv')

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
# cell type distribution 
get_cell<-function(working_dataset, sig_table, origin_gene, origin_tissue,col_scheme){
  #read in sc-seq matrix 
  sc_matrix = sc_matrix
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
  top_genes2 = cor_table_sc[cor_table_sc$pvalue<0.01,]
  top_genes2$color = col_scheme[match(top_genes2$tissue, names(col_scheme))]
  top_genes2 = top_genes2[order(abs(top_genes2$bicor), decreasing = T),]
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