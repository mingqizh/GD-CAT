
##################################################
#SQL tables to start
library(stats)
library(WGCNA)
library(ggpubr)
library(reshape2)
allowWGCNAThreads()
load('mouse.RData')

######################################
#Mingqi These scripts are just reproducing the adipoQ file that you have so ignore as these will be SQL pulls
working_dataset<-HF
tissue2 <- working_dataset

tissue1 <- working_dataset[, colnames(working_dataset) == 'Adipoq_adipose']

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
origin_gene = 'Adipoq'
origin_tissue = 'adipose'
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


#I set the tissue manually but here is where the user would select the tissue based on what you showed me in clicking the pie chart
select_tissue = 'adipose'
pie_bin = 0.1

setEnrichrSite("Enrichr")
dbs1 <- c("GO_Biological_Process_2021", "MGI_Mammalian_Phenotype_Level_4_2019", "Reactome_2022")

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

plotEnrich(enriched[[2]], showTerms = 8, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('positive gene correlations with ', origin_gene_tissue, ' ', names(enriched[2])))

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

plotEnrich(enriched[[2]], showTerms = 8, numChar = 40, y = "Count", orderBy = "P.value") + ggtitle(paste0('negative gene correlations with ', origin_gene_tissue, ' ', names(enriched[2])))

plotEnrich(enriched[[3]], showTerms = 4, numChar = 45, y = "Count", orderBy = "P.value") + ggtitle(paste0('negative gene correlations with ', origin_gene_tissue, ' ', names(enriched[3])))


#Also insert other rank-based pathways
library(enrichR)
#Note this will be library(org.Mm.eg.db) for mouse
library(org.Hs.eg.db)
library("pathview")
library(clusterProfiler)
library(enrichplot)
library(cowplot)
binned_sig_prots= sig_table %>%
  dplyr::group_by(qcat, tissue_2) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))

tissue_freqs = binned_sig_prots %>% dplyr::group_by(qcat) %>% dplyr::summarise(sum(n))

binned_sig_prots$tot_count = tissue_freqs$`sum(n)`[match(binned_sig_prots$qcat, tissue_freqs$qcat)]
binned_sig_prots$qcat1 = paste0(binned_sig_prots$qcat, ' ', binned_sig_prots$tot_count, ' genes')

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

gse@result$Description[1:20]
##emapplot pathway networks

x2 <- pairwise_termsim(gse)
emapplot(x2, showCategory = 20)+ ggtitle("Relationship between the top 20 most significantly GSE - GO terms (padj.)")
goplot(gse)




#####Clinical trait cors
full_cors1 = bicorAndPvalue(tissue1, working_traits, use = 'p')
cor_table1 = reshape2::melt(full_cors1$bicor)
new_p = reshape2::melt(full_cors1$p)
cor_table1$Var1=NULL
colnames(cor_table1) = c( 'trait_name', 'bicor')
#can drop here to clear CPU
full_cors1=NULL

cor_table1$pvalue = signif(new_p$value, 3)
cor_table1 = cor_table1[order(cor_table1$pvalue, decreasing = F),]

#we should let users export this table!



#have user select the top traits to visualize
topn=30
cc1 = cor_table1[1:topn,]
cc1$direction = ifelse(cc1$bicor>0, 'positive correlation', 'negative correlation')
cc1$logp = -log10(cc1$pvalue)
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
) + ylab('-log10(pvalue)') + xlab('') + ggtitle(paste0(origin_gene_tissue, ' ~ trait correlations'))

#have user enter the gene ~ trait correlation.  Similar to gene ~ gene correlations


####################################################
#next we will generate the network.  I run the analyses but the SQL pull should be based on this vector:


number_orig_gene = 20
number_peripheral_genes = 20


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
table(map1$cols)
map2 = reshape2::dcast(map1, gene_tissue_1 ~ gene_tissue_2, value.var = 'bicor', fun.aggregate = mean)
row.names(map2) = map2$gene_tissue_1
map2$gene_tissue_1 = NULL
cols_set = as.data.frame(row.names(map2))
row.names(cols_set) = row.names(map2)
cols_set$cols = map1$cols[match(row.names(cols_set), map1$gene_tissue_1)]
cols_set$cols = ifelse(row.names(cols_set) %in% origin_gene_tissue, 'gray5', paste0(cols_set$cols))
head(map2)
qgraph(map2, minimum = 0.4, cut = 0.9, vsize = 3, color=cols_set$cols, labels=gsub("\\_.*","",row.names(map2)), legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=2, directed=F, labels = colnames(map2))

qgraph(map2, minimum = 0.3, cut = 0.8, vsize = 3, color=cols_set$cols, legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3",  directed=F, labels = F) + ggtitle(paste0('Undirected network for ', origin_gene_tissue, ' without labels')) 

dev.off()


