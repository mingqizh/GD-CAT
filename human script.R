##################################################
#For Natalie
#produce the initial tables to upload for SQL
library(stats)
library(WGCNA)
library(reshape2)
allowWGCNAThreads()
setwd('G:/My Drive/lab files/endocrine signaling app')
load('G:/My Drive/lab files/sex-difference myokine study/GTEx NA included env.RData')
GTEx_full=NULL
working_dataset=female
row.names(working_dataset) = working_dataset$gene_tissue
working_dataset$gene_tissue=NULL
working_dataset = as.data.frame(t(working_dataset))
GTEx_subfiltered=NULL



######################################
#Mingqi These scripts are just reproducing the adipoQ file that you have so ignore as these will be SQL pulls
tissue2 <- working_dataset[,grepl('Adipose - Subcutaneous', colnames(working_dataset)) | grepl('Adipose - Visceral (Omentum)', colnames(working_dataset), fixed=T) | grepl('Brain - Hypothalamus', colnames(working_dataset)) | grepl('Brain - Hippocampus', colnames(working_dataset))  | grepl('Small Intestine - Terminal Ileum', colnames(working_dataset), fixed=T) | grepl('Stomach', colnames(working_dataset), fixed=T) | grepl('Thyroid', colnames(working_dataset), fixed=T) | grepl('Pancreas', colnames(working_dataset), fixed=T) | grepl('Spleen', colnames(working_dataset), fixed=T) | grepl('Muscle - Skeletal', colnames(working_dataset), fixed=T) | grepl('Pituitary', colnames(working_dataset), fixed=T) | grepl('Artery - Coronary', colnames(working_dataset), fixed=T) | grepl('Liver', colnames(working_dataset), fixed=T) | grepl('Kidney - Cortex', colnames(working_dataset), fixed=T) | grepl('Heart - Left Ventricle', colnames(working_dataset), fixed=T) | grepl('Colon - Transverse', colnames(working_dataset), fixed=T) | grepl('Colon - Sigmoid', colnames(working_dataset), fixed=T) | grepl('Adrenal Gland', colnames(working_dataset), fixed=T) |  grepl('Artery - Aorta', colnames(working_dataset), fixed=T),]


tissue1 <- working_dataset[, colnames(working_dataset) == 'ARNTL_Colon - Transverse']


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
origin_gene = 'ARNTL'
origin_tissue = 'Colon - Transverse'
origin_gene_tissue = paste0(origin_gene, '_', origin_tissue)

sig_table = sql_pull_1[sql_pull_1$qvalue<0.1,]
sig_table$qcat =ifelse(sig_table$qvalue<0.01, 'q<0.01', 'q<0.1')
sig_table$qcat =ifelse(sig_table$qvalue<0.0001, 'q<0.0001', paste0(sig_table$qcat))

##set the color scheme up front
col_scheme = rev(met.brewer('Austria', length(unique(sig_table$tissue_2))))
names(col_scheme) = unique(sig_table$tissue_2)

#Up front I think it will be helpful to have the user see the top-ranked genes correlating both within and excluding origin.  I set 30 as max_gene_length but this might be cool to let user enter.  We should max out at 40 or 50
max_gene_length = 40

top_genes = sig_table[1:max_gene_length,]
top_genes$color = col_scheme[match(top_genes$tissue_2, names(col_scheme))]


ggplot(top_genes, aes(x=fct_reorder(gene_tissue_2, abs(bicor), .desc=T), y=bicor, fill=tissue_2)) + geom_col(position = position_dodge2(reverse = TRUE)) +  theme_classic() + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + scale_fill_manual(values = top_genes$color[order(abs(top_genes$bicor), decreasing = T)]) + xlab('') + ylab('bicor coefficent') + ggtitle(paste0('Top genes correlated with ', origin_gene_tissue))


#read in sc-seq matrix 
sc_matrix = read.csv('G:/My Drive/lab files/endocrine signaling app/sc_deconvolution/full pan-tissue DECON abundances.csv')

#filter sc-seq matrix based on working_dataset object 
row.names(sc_matrix) = sc_matrix$GTExID
sc_matrix$GTExID = NULL
sc_matrix = sc_matrix[row.names(sc_matrix) %in% row.names(working_dataset),]
tissue1 <- working_dataset[,grepl('ARNTL_Colon - Transverse', colnames(working_dataset), fixed=T)]

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
cor_table_sc$cell_type = gsub("\\_.*","",cor_table_sc$cell_tissue)
head(cor_table_sc)
top_genes2 = cor_table_sc[cor_table_sc$tissue==origin_tissue,]
top_genes2$color = col_scheme[match(top_genes2$tissue, names(col_scheme))]
top_genes2 = top_genes2[order(abs(top_genes2$bicor), decreasing = T),]

#pdf(file = 'significant cell types aadiponectin cors.pdf')
ggplot(top_genes2, aes(x=fct_reorder(cell_type, abs(bicor), .desc=T), y=bicor, fill=color)) + geom_col(position = position_dodge2(reverse = TRUE)) +  theme_classic() + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + scale_fill_manual(values = top_genes2$color[order(abs(top_genes2$bicor), decreasing = T)]) + xlab('') + ylab('bicor coefficent') + ggtitle(paste0('Significant cell-types correlated with ', origin_gene_tissue, 'P < 1e-3'))
dev.off()



top_genes2 = cor_table_sc[!cor_table_sc$tissue==origin_tissue,]
top_genes2 = top_genes2[!top_genes2$tissue=='Lung',]

top_genes2$color = col_scheme[match(top_genes2$tissue, names(col_scheme))]
top_genes2 = top_genes2[order(abs(top_genes2$bicor), decreasing = T),]
top_genes2 = top_genes2[1:20,]

#pdf(file = 'significant cell types aadiponectin cors.pdf')
ggplot(top_genes2, aes(x=fct_reorder(cell_tissue, abs(bicor), .desc=T), y=bicor, fill=color)) + geom_col(position = position_dodge2(), fill =fct_reorder(top_genes2$color, abs(top_genes2$bicor), .desc=T) ) +  theme_classic() + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+ xlab('') + ylab('bicor coefficent') + ggtitle(paste0('Significant cell-types correlated with ', origin_gene_tissue, 'P < 1e-3'))
dev.off()


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
  select_tissue = 'Adipose - Subcutaneous'
pie_bin = 0.1


library(enrichR)
#Note this will be library(org.Mm.eg.db) for mouse
library(org.Hs.eg.db)
library("pathview")
library(clusterProfiler)
library(enrichplot)
library(cowplot)

tissue_table = binned_sig_prots[binned_sig_prots$qcat =='q<0.1',]
select_tissue = tissue_table$tissue_2[1]
pie_bin = 0.1
pp1 = sig_table[sig_table$tissue_2 %in% select_tissue,]
pp1 = pp1[pp1$pvalue < pie_bin,]
head(pp1)
fc_dko = scales::rescale(pp1$bicor, to=c(-3, 3))

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
  pvalueCutoff = 0.5,
  pAdjustMethod = "BH") 
str(gse)

## gsego https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/gseGO
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign) +
  ggtitle(paste0('GSEA pathways from positive and negative correlations ',  origin_gene_tissue, ' in ', select_tissue))

##emapplot pathway networks

x2 <- pairwise_termsim(gse)
emapplot(x2, showCategory = 20)+ ggtitle("Relationship between the top 20 most significantly GSE - GO terms (padj.)")
goplot(gse)



#we can remove these

setEnrichrSite("Enrichr")
dbs1 <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2022")

#first for positive
pp1 = sig_table[sig_table$tissue_2 %in% select_tissue,]
head(pp1)
pp1 = pp1[pp1$qvalue < pie_bin,]
pp1 = pp1[pp1$bicor>0,]
pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
pp2 = pp1[1:pp1_length,]
gg1 = pp2$gene_symbol_2
enriched <- enrichr(gg1, dbs1)
plotEnrich(enriched[[1]], showTerms = 4, numChar = 20, y = "Count", orderBy = "P.value") + ggtitle(paste0('positive gene correlations with ', origin_gene_tissue, ' ', names(enriched[1])))

#plotEnrich(enriched[[2]], showTerms = 8, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('positive gene correlations with ', origin_gene_tissue, ' ', names(enriched[2])))

plotEnrich(enriched[[3]], showTerms = 4, numChar = 45, y = "Count", orderBy = "P.value") + ggtitle(paste0('positive gene correlations with ', origin_gene_tissue, ' ', names(enriched[3])))


#now for negative
pp1 = sig_table[sig_table$tissue_2 %in% select_tissue,]
pp1 = pp1[pp1$bicor<0,]
pp1 = pp1[pp1$qvalue < pie_bin,]
pp1_length = ifelse(length(row.names(pp1)) > 500, as.numeric(500), as.numeric(length(row.names(pp1))))
pp2 = pp1[1:pp1_length,]
gg1 = pp2$gene_symbol_2
enriched <- enrichr(gg1, dbs1)
plotEnrich(enriched[[1]], showTerms = 4, numChar = 20, y = "Count", orderBy = "P.value") + ggtitle(paste0('negative gene correlations with ', origin_gene_tissue, ' ', names(enriched[1])))

#plotEnrich(enriched[[2]], showTerms = 8, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('negative gene correlations with ', origin_gene_tissue, ' ', names(enriched[2])))

plotEnrich(enriched[[3]], showTerms = 4, numChar = 45, y = "Count", orderBy = "P.value") + ggtitle(paste0('negative gene correlations with ', origin_gene_tissue, ' ', names(enriched[3])))






####################################################
#next we will generate the network.  I run the analyses but the SQL pull should be based on this vector:
paths_table = read.delim('G:/My Drive/Datasets/Human/genome files/uniprot-human-genes and goterms mapping.tab')
head(sig_table)
select_paths = paths_table[grepl('oxidation|fatty acid', paths_table$Gene.ontology..biological.process.),]
select_paths1 = paths_table[grepl('contraction|sarcomer|muscle', paths_table$Gene.ontology..biological.process.),]

number_orig_gene = 20
number_peripheral_genes = 20


origin_pull = sig_table[sig_table$tissue_2==origin_tissue,]
orig_network_genes = as.vector(origin_pull$gene_tissue_2[1:number_orig_gene])
#peripheral_pull = sig_table[!sig_table$tissue_2==origin_tissue,]
peripheral_pull = sig_table[sig_table$tissue_2=='Muscle - Skeletal',]
periph_network_genes = as.vector(peripheral_pull$gene_tissue_2[1:number_peripheral_genes])

sql_pull_list = c(origin_gene_tissue, orig_network_genes, periph_network_genes)
#sql_pull_list = c(origin_gene_tissue, orig_network_genes)
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
table(map1$cols)
map1$cols = ifelse(map1$gene_symbol_2 %in% select_paths$Gene.names...primary.., 'firebrick2', paste0(map1$cols))
map1$cols = ifelse(map1$gene_symbol_2 %in% select_paths1$Gene.names...primary.., 'darkorchid1', paste0(map1$cols))
map2 = reshape2::dcast(map1, gene_tissue_1 ~ gene_tissue_2, value.var = 'bicor', fun.aggregate = mean)
row.names(map2) = map2$gene_tissue_1
map2$gene_tissue_1 = NULL
cols_set = as.data.frame(row.names(map2))
row.names(cols_set) = row.names(map2)
cols_set$cols = map1$cols[match(row.names(cols_set), map1$gene_tissue_1)]
cols_set$cols = ifelse(row.names(cols_set) %in% origin_gene_tissue, 'gray5', paste0(cols_set$cols))
head(map2)
qgraph(map2, minimum = 0.4, cut = 0.9, vsize = 3, color=cols_set$cols, labels=gsub("\\_.*","",row.names(map2)), legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=2, directed=F, labels = colnames(map2)) + ggtitle(paste0('Undirected network for ', origin_gene_tissue, ' with labels')) 

qgraph(map2, minimum = 0.3, cut = 0.8, vsize = 3, color=cols_set$cols, legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3",  directed=F, labels = F) + ggtitle(paste0('Undirected network for ', origin_gene_tissue, ' without labels')) 



###############################Age Seq BMI
origin_gene = 'NOS3'
origin_tissue = 'Artery - Aorta'
origin_gene_tissue = paste0(origin_gene, '_', origin_tissue)
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


library(MetBrewer)
#######################gene ~ gene cor

select_gene2 = 'ERFE_Muscle - Skeletal'
scat_gene_set1 = working_dataset[,  colnames(working_dataset) %in% origin_gene_tissue]
scat_gene_set2 = working_dataset[,  colnames(working_dataset) %in% select_gene2]

scat_df = as.data.frame(cbind(scat_gene_set1, scat_gene_set2))
colnames(scat_df) = c('sgene_1', 'sgene2')
cc1 = bicorAndPvalue(scat_df$sgene2, scat_df$sgene_1)
my_palette <- met.brewer("Peru1")
ggplot(scat_df, aes(x=sgene_1, y=sgene2)) + theme_classic() +  geom_hex(bins = 100) + scale_fill_distiller(palette = "Blues", direction=-1)+ geom_smooth(method = 'lm', color = "darkorange")+  ggtitle(paste0(origin_gene_tissue, ' ~ ', select_gene2, ' r=', round(cc1$bicor, 3), ' p=', signif(cc1$p, 3))) + xlab(paste0('Normalized ', origin_gene_tissue, ' expression')) + ylab(paste0('Normalized ', select_gene2, ' expression')) 
dev.off()
