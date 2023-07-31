library(shiny)
library(shinyWidgets)
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
library(data.table)
library(networkD3)
##################################################
#For Natalie
#produce the initial tables to upload for SQL
library(stats)
library(WGCNA)
library(reshape2)
allowWGCNAThreads()
setwd('G:/My Drive/lab files/endocrine signaling app')
load("C:/Users/mingqiz7/Desktop/GTEx app/data/GTEx env for endocrine app - Both sexes combined.RData")
GTEx_full=NULL
working_dataset=GTEx_subfiltered
row.names(working_dataset) = working_dataset$gene_tissue
working_dataset$gene_tissue=NULL
working_dataset = as.data.frame(t(working_dataset))
working_dataset[1:5,1:5]
GTEx_subfiltered=NULL


#These list all the tissues available that will be correlated in a pairwise fashion.  
#tissue1 <- working_dataset[,grepl('Adipose - Subcutaneous', colnames(working_dataset)) | grepl('Adipose - Visceral (Omentum)', colnames(working_dataset), fixed=T) | grepl('Brain - Hypothalamus', colnames(working_dataset)) | grepl('Brain - Hippocampus', colnames(working_dataset)) | grepl('Lung', colnames(working_dataset)) | grepl('Small Intestine - Terminal Ileum', colnames(working_dataset), fixed=T) | grepl('Stomach', colnames(working_dataset), fixed=T) | grepl('Thyroid', colnames(working_dataset), fixed=T) | grepl('Pancreas', colnames(working_dataset), fixed=T) | grepl('Spleen', colnames(working_dataset), fixed=T) | grepl('Muscle - Skeletal', colnames(working_dataset), fixed=T) | grepl('Pituitary', colnames(working_dataset), fixed=T) | grepl('Artery - Coronary', colnames(working_dataset), fixed=T) | grepl('Liver', colnames(working_dataset), fixed=T) | grepl('Kidney - Cortex', colnames(working_dataset), fixed=T) | grepl('Heart - Left Ventricle', colnames(working_dataset), fixed=T) | grepl('Colon - Transverse', colnames(working_dataset), fixed=T) | grepl('Colon - Sigmoid', colnames(working_dataset), fixed=T) | grepl('Adrenal Gland', colnames(working_dataset), fixed=T) |  grepl('Artery - Aorta', colnames(working_dataset), fixed=T),]


#Here is an example for one tissue with only 100 genes from another as a toy example (note that tissues will be correlated with themselves as well)
tissue1 <- working_dataset[, grepl('Adipose - Subcutaneous', colnames(working_dataset))  ]
tissue1 = na.omit(tissue1)
tissue2 <- working_dataset[, grepl('Adipose - Visceral (Omentum)', colnames(working_dataset), fixed=T) ]

tissue2 = tissue2[,1:1000]
tissue2 = na.omit(tissue2)

#make sure the individuals are matching
tissue1 = tissue1[row.names(tissue1) %in% row.names(tissue2),]
tissue2 = tissue2[row.names(tissue2) %in% row.names(tissue1),]

full_cors = bicorAndPvalue(tissue1, tissue2, use = 'p')
cor_table = reshape2::melt(full_cors$bicor)
new_p = reshape2::melt(full_cors$p)

colnames(cor_table) = c('gene_tissue_1', 'gene_tissue_2', 'bicor')
#can drop here to clear CPU
full_cors=NULL

cor_table$pvalue = signif(new_p$value, 3)
cor_table$bicor = round(cor_table$bicor, 3)
cor_table$qvalue = signif(p.adjust(cor_table$pvalue, "BH"), 3)
cor_table = cor_table[order(cor_table$qvalue, decreasing=F),]
cor_table = na.omit(cor_table)
cor_table$gene_symbol_1 = gsub("\\_.*","",cor_table$gene_tissue_1)
cor_table$tissue_1 = gsub(".*_","",cor_table$gene_tissue_1)
cor_table = cor_table[!is.na(cor_table$tissue_1),]

cor_table$gene_symbol_2 = gsub("\\_.*","",cor_table$gene_tissue_2)
cor_table$tissue_2 = gsub(".*_","",cor_table$gene_tissue_2)
cor_table = cor_table[!is.na(cor_table$tissue_2),]

#I wrote as a csv but lets check with mingqi to see which format is more amenable
write.csv(cor_table, file = 'two tissue pairs example.csv', row.names=F)



######################################
#Mingqi These scripts are just reproducing the adipoQ file that you have so ignore as these will be SQL pulls
tissue2 <- working_dataset[,grepl('Adipose - Subcutaneous', colnames(working_dataset)) | grepl('Adipose - Visceral (Omentum)', colnames(working_dataset), fixed=T) | grepl('Brain - Hypothalamus', colnames(working_dataset)) | grepl('Brain - Hippocampus', colnames(working_dataset)) | grepl('Lung', colnames(working_dataset)) | grepl('Small Intestine - Terminal Ileum', colnames(working_dataset), fixed=T) | grepl('Stomach', colnames(working_dataset), fixed=T) | grepl('Thyroid', colnames(working_dataset), fixed=T) | grepl('Pancreas', colnames(working_dataset), fixed=T) | grepl('Spleen', colnames(working_dataset), fixed=T) | grepl('Muscle - Skeletal', colnames(working_dataset), fixed=T) | grepl('Pituitary', colnames(working_dataset), fixed=T) | grepl('Artery - Coronary', colnames(working_dataset), fixed=T) | grepl('Liver', colnames(working_dataset), fixed=T) | grepl('Kidney - Cortex', colnames(working_dataset), fixed=T) | grepl('Heart - Left Ventricle', colnames(working_dataset), fixed=T) | grepl('Colon - Transverse', colnames(working_dataset), fixed=T) | grepl('Colon - Sigmoid', colnames(working_dataset), fixed=T) | grepl('Adrenal Gland', colnames(working_dataset), fixed=T) |  grepl('Artery - Aorta', colnames(working_dataset), fixed=T),]


tissue1 <- working_dataset[,grepl('ADIPOQ_Adipose - Subcutaneous', colnames(working_dataset), fixed=T)]


#tissue1 = tissue1[row.names(tissue1) %in% row.names(tissue2),]
#tissue2 = tissue2[row.names(tissue2) %in% row.names(tissue1),]

full_cors = bicorAndPvalue(tissue1, tissue2, use = 'p')
cor_table = reshape2::melt(full_cors$bicor)
new_p = reshape2::melt(full_cors$p)

colnames(cor_table) = c('gene_tissue_1', 'gene_tissue_2', 'bicor')
#can drop here to clear CPU
full_cors=NULL

cor_table$pvalue = signif(new_p$value, 3)
cor_table$bicor = round(cor_table$bicor, 3)
cor_table$qvalue = signif(p.adjust(cor_table$pvalue, "BH"), 3)
cor_table = cor_table[order(cor_table$qvalue, decreasing=F),]
cor_table = na.omit(cor_table)
cor_table$gene_symbol_1 = gsub("\\_.*","",cor_table$gene_tissue_1)
cor_table$tissue_1 = gsub(".*_","",cor_table$gene_tissue_1)
cor_table = cor_table[!is.na(cor_table$tissue_1),]

cor_table$gene_symbol_2 = gsub("\\_.*","",cor_table$gene_tissue_2)
cor_table$tissue_2 = gsub(".*_","",cor_table$gene_tissue_2)
cor_table = cor_table[!is.na(cor_table$tissue_2),]

sql_pull_1 = cor_table





#Mingqi here is where Ill start for the app.  I called the initial object that you will retrieve from SQL sql_pull_1.  These apps should be all that you need also
library(ggplot2)
library(qgraph)
library(dplyr)
library(forcats)
library(enrichR)
library(MetBrewer)
#here I set the origin gene_tissue manually, but the user will enter 
origin_gene = 'ADIPOQ'
origin_tissue = 'Adipose - Subcutaneous'
origin_gene_tissue = paste0(origin_gene, '_', origin_tissue)

a<-fread('significant crosstissue enrichments ADIPOQ_Adipose - Subcutaneous - Final version.csv')
a<-cor_table
a$gene_tissue_1 = paste0(a$gene_symbol_1, '_', a$tissue_1)
a$gene_tissue_2 = paste0(a$gene_symbol_2, '_', a$tissue_2)
sig_table = a[a$qvalue<0.1,]
sig_table$qcat =ifelse(sig_table$qvalue<0.01, 'q<0.01', 'q<0.1')
sig_table$qcat =ifelse(sig_table$qvalue<0.0001, 'q<0.0001', paste0(sig_table$qcat))

##set the color scheme up front
col_scheme = rev(met.brewer('Austria', length(unique(sig_table$tissue_2))))
names(col_scheme) = unique(sig_table$tissue_2)

## sex age
library(gridExtra)
library(ggpubr)
met_dat = read.delim('GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')

sex_table = read.delim('GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
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

#Up front I think it will be helpful to have the user see the top-ranked genes correlating both within and excluding origin.  I set 30 as max_gene_length but this might be cool to let user enter.  We should max out at 40 or 50
max_gene_length = 40
top_genes = sig_table[1:max_gene_length,]
top_genes$color = col_scheme[match(top_genes$tissue_2, names(col_scheme))]

ggplot(top_genes, aes(x=fct_reorder(gene_tissue_2, abs(bicor), .desc=T), y=bicor, fill=tissue_2)) + geom_col(position = position_dodge2(reverse = TRUE)) +  theme_classic() + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values = top_genes$color[order(abs(top_genes$bicor), decreasing = T)]) + xlab('') + ylab('bicor coefficent') + ggtitle(paste0('Top genes correlated with ', origin_gene_tissue))

#same but remove origin
sig_table1 = sig_table[!sig_table$tissue_2 %in% origin_tissue,]
top_genes = sig_table1[1:max_gene_length,]
top_genes$color = col_scheme[match(top_genes$tissue_2, names(col_scheme))]

ggplot(top_genes, aes(x=fct_reorder(gene_tissue_2, abs(bicor), .desc=T), y=bicor, fill=tissue_2)) + geom_col(position = position_dodge2(reverse = TRUE)) +  theme_classic() + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values = top_genes$color[order(abs(top_genes$bicor), decreasing = T)]) + xlab('') + ylab('bicor coefficent') + ggtitle(paste0('Top genes correlated with ', origin_gene_tissue, ' origin tissue removed'))


#produce pie charts including origin tissue
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

ggplot(binned_sig_prots, aes(x = "", y = freq, fill =tissue_2)) + 
  geom_bar(stat = "identity", width = 1, position = position_fill()) + blank_theme + theme(plot.title=element_text(size=25, face="bold")) +
  theme(axis.text.x=element_blank())+ scale_fill_manual(values=col_scheme[match(binned_sig_prots$tissue_2, names(col_scheme))]) +
  coord_polar(theta = "y") + 
  facet_wrap( ~ qcat1)


#############################################
#produce the pie charts excluding origin tissue
sig_table1 = sig_table[!sig_table$tissue_2 %in% origin_tissue,]

binned_sig_prots= sig_table1 %>%
  dplyr::group_by(qcat, tissue_2) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))

tissue_freqs = binned_sig_prots %>% dplyr::group_by(qcat) %>% dplyr::summarise(sum(n))

binned_sig_prots$tot_count = tissue_freqs$`sum(n)`[match(binned_sig_prots$qcat, tissue_freqs$qcat)]
binned_sig_prots$qcat1 = paste0(binned_sig_prots$qcat, ' ', binned_sig_prots$tot_count, ' genes')

ggplot(binned_sig_prots, aes(x = "", y = freq, fill =tissue_2)) + 
  geom_bar(stat = "identity", width = 1, position = position_fill()) + blank_theme + theme(plot.title=element_text(size=25, face="bold")) +
  theme(axis.text.x=element_blank())+ scale_fill_manual(values=col_scheme[match(binned_sig_prots$tissue_2, names(col_scheme))]) +
  coord_polar(theta = "y") + 
  facet_wrap( ~ qcat1)




################################################
#run the pathways corresponding to each tissue.  I revised these to bin separately by negatively vs positively correlated which I think will be important for the user

#I set the tissue manually but here is where the user would select the tissue based on what you showed me in clicking the pie chart
select_tissue = 'Muscle - Skeletal'
pie_bin = 0.01

setEnrichrSite("Enrichr")
dbs1 <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2022")

#first for positive
pp1 = sig_table[sig_table$tissue_2 %in% select_tissue,]
pp1 = pp1[pp1$qvalue < pie_bin,]
pp1 = pp1[pp1$bicor>0,]
pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
pp2 = pp1[1:pp1_length,]
gg1 = pp2$gene_symbol_2
enriched <- enrichr(gg1, dbs1)
plotEnrich(enriched[[1]], showTerms = 8, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('positive gene correlations with ', origin_gene_tissue, ' ', names(enriched[1])))

plotEnrich(enriched[[2]], showTerms = 8, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('positive gene correlations with ', origin_gene_tissue, ' ', names(enriched[2])))

plotEnrich(enriched[[3]], showTerms = 8, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('positive gene correlations with ', origin_gene_tissue, ' ', names(enriched[3])))


#now for negative
pp1 = sig_table[sig_table$tissue_2 %in% select_tissue,]
pp1 = pp1[pp1$bicor<0,]
pp1 = pp1[pp1$qvalue < pie_bin,]
pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
pp2 = pp1[1:pp1_length,]
gg1 = pp2$gene_symbol_2
enriched <- enrichr(gg1, dbs1)
plotEnrich(enriched[[1]], showTerms = 8, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('negative gene correlations with ', origin_gene_tissue, ' ', names(enriched[1])))

plotEnrich(enriched[[2]], showTerms = 8, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('negative gene correlations with ', origin_gene_tissue, ' ', names(enriched[2])))

plotEnrich(enriched[[3]], showTerms = 8, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('negative gene correlations with ', origin_gene_tissue, ' ', names(enriched[3])))




####################################################
#next we will generate the network.  I run the analyses but the SQL pull should be based on this vector:
network_genes = sig_table %>%
  group_by(tissue_2) %>%
  top_n(200, gene_tissue_2)

sql_pull_list = c(as.vector(network_genes$gene_tissue_2), origin_gene_tissue)

######################
#I reran cors for these but skip to line 283
tissue1 <- working_dataset[, colnames(working_dataset)  %in% sql_pull_list]
full_cors = bicorAndPvalue(tissue1, tissue1, use = 'p')
cor_table = reshape2::melt(full_cors$bicor)
new_p = reshape2::melt(full_cors$p)
colnames(cor_table) = c('gene_tissue_1', 'gene_tissue_2', 'bicor')
#can drop here to clear CPU
full_cors=NULL
cor_table$pvalue = signif(new_p$value, 3)
cor_table$bicor = round(cor_table$bicor, 3)
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

######################################
#the new SQL pull table is not listed as sql_pull_2
#this will take work on your end to interface with pulling the correct genes from tissue_1 based on the number selected in the app.  I arbitrarily chose 10 here and you obviously wouldn't jsut apply top_n since it will differ.  The final list is called network_gene_lsit
library(data.table)
sql_pull_2<-fread(file = 'sql_pull_2.csv')
network_gene_list = sql_pull_2 %>%  
  group_by(tissue_1) %>%
  slice(1:10)


network_plot_table = sql_pull_2[sql_pull_2$gene_tissue_1 %in% network_gene_list$gene_tissue_1,]
network_plot_table = network_plot_table[network_plot_table$gene_tissue_2 %in% network_gene_list$gene_tissue_1,]

map1 = network_plot_table
map1$cols = col_scheme[match(map1$tissue_1, names(col_scheme))]
map2 = dcast(map1, gene_tissue_1 ~ gene_tissue_2, value.var = 'bicor', fun.aggregate = mean)
row.names(map2) = map2$gene_tissue_1
map2$gene_tissue_1 = NULL
cols_set = as.data.frame(row.names(map2))
row.names(cols_set) = row.names(map2)
cols_set$cols = map1$cols[match(row.names(cols_set), map1$gene_tissue_1)]
cols_set$cols = ifelse(row.names(cols_set) %in% origin_gene_tissue, 'gray5', paste0(cols_set$cols))

qgraph(map2, minimum = 0.3, cut = 0.8, vsize = 3, color=cols_set$cols, legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=3, directed=F, labels = colnames(map2)) + ggtitle(paste0('Undirected network for ', origin_gene_tissue, ' with labels')) 

qgraph(map2, minimum = 0.3, cut = 0.8, vsize = 4, color=cols_set$cols, legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3",  directed=F, labels = F) + ggtitle(paste0('Undirected network for ', origin_gene_tissue, ' without labels')) 

number_orig_gene = 50
number_peripheral_genes = 150
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

qgraph(map2, minimum = 0.3, cut = 0.8, vsize = 3, color=cols_set$cols, legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=6, directed=F, labels = colnames(map2)) 
qgraph(map2, minimum = 0.3, cut = 0.8, vsize = 3, color=cols_set$cols, legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3",  directed=F, labels = F) + ggtitle(paste0('Undirected network for ', origin_gene_tissue, ' without labels')) 
