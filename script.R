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
working_dataset=fread("first 100 gene cocorrelation across two adipose tissues - melted data for SQL input.csv")
working_dataset<-annot
working_dataset$gene_tissue_1 = paste0(working_dataset$gene_symbol_1, '_', working_dataset$tissue_1)
working_dataset$gene_tissue_2 = paste0(working_dataset$gene_symbol_2, '_', working_dataset$tissue_2)
a<- working_dataset$gene_tissue_2
working_dataset$gene_tissue_2=NULL
working_dataset = as.data.frame(t(working_dataset))
colnames(working_dataset)<-a
test1 = working_dataset[,grepl('AARS2', colnames(working_dataset))]
colnames(test1)

gene_1 = 'AARS2'
origin_tissue = 'Adipose - Subcutaneous'
gene_tissue1 = paste0(gene_1, '_', origin_tissue)
tissue1 <- working_dataset[,grepl('Adipose - Subcutaneous', colnames(working_dataset))
                           | grepl('Adipose - Visceral (Omentum)', colnames(working_dataset), fixed=T)
                           | grepl('Brain - Hypothalamus', colnames(working_dataset), fixed=T)
                           | grepl('Colon - Transverse', colnames(working_dataset), fixed=T)
                           | grepl('Spleen', colnames(working_dataset), fixed=T)
                           | grepl('Small Intestine - Terminal Ileum', colnames(working_dataset), fixed=T)
                           | grepl('Artery - Coronary', colnames(working_dataset), fixed=T)
                           | grepl('Stomach', colnames(working_dataset), fixed=T)
                           | grepl('Thyroid', colnames(working_dataset), fixed=T)
                           | grepl('Pancreas', colnames(working_dataset), fixed=T)
                           | grepl('Muscle - Skeletal', colnames(working_dataset), fixed=T)
                           | grepl('Pituitary', colnames(working_dataset), fixed=T)
                           | grepl('Liver', colnames(working_dataset), fixed=T)
                           | grepl('Kidney - Cortex', colnames(working_dataset), fixed=T)
                           | grepl('Heart - Left Ventricle', colnames(working_dataset), fixed=T)
                           | grepl('Colon - Sigmoid', colnames(working_dataset), fixed=T)
                           | grepl('Adrenal Gland', colnames(working_dataset), fixed=T) 
                           | grepl('Artery - Aorta', colnames(working_dataset), fixed=T)
                           | grepl('Brain - Hippocampus', colnames(working_dataset), fixed=T)
                           | grepl('Lung', colnames(working_dataset), fixed=T)]
tissue1 = as.data.frame(tissue1)

origin = tissue1[,grepl(gene_tissue1, colnames(tissue1))]
target = tissue1[,!grepl(gene_tissue1, colnames(tissue1))]

full_cors = bicorAndPvalue(origin, target, use = 'p')
cor_table = reshape2::melt(full_cors$bicor)
cor_table$Var1=NULL
colnames(cor_table) = c('gene_tissue', 'bicor')
new_p = reshape2::melt(full_cors$p)

cor_table$pvalue = new_p$value[match(cor_table$gene_tissue, new_p$Var2)]
cor_table = na.omit(cor_table)
qest = qvalue(cor_table$pvalue, pi0 = 1)
cor_table$qvalue = qest$qvalues
cor_table$gene_symbol = gsub("\\_.*","",cor_table$gene_tissue)
cor_table$tissue = gsub(".*_","",cor_table$gene_tissue)
cor_table = cor_table[!is.na(cor_table$tissue),]
cor_table = cor_table[order(cor_table$pvalue, decreasing = F),]

res1 = cor_table[cor_table$pvalue<0.01,]
res1 = na.omit(res1)
write.csv(res1, file = paste0('significant crosstissue enrichments ', gene_tissue1, '.csv'), row.names = F)

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

binned_sig_prots$tissue_2 <- factor(binned_sig_prots$tissue_2, levels = levels(binned_sig_prots$n[binned_sig_prots$qcat=="q<0.01 11557 genes",]))

pdf(file = paste0('crosstissue gene enrichments ', gene_tissue1, '.pdf'))
ggplot(binned_sig_prots, aes(x = "", y = freq, fill =tissue_2)) + 
  geom_bar(stat = "identity", width = 1, position = position_fill()) + blank_theme + theme(plot.title=element_text(size=25, face="bold")) +
  theme(axis.text.x=element_blank())+ scale_fill_manual(values=col_scheme) +
  coord_polar(theta = "y") + 
  facet_wrap( ~ qcat1)
dev.off()


## pathway tissue2
working_dataset=annot
tissue_set = unique(gsub(".*_","",annot$tissue_2))
new_cols = met.brewer('Moreau', length(tissue_set))
names(new_cols) = unique(tissue_set)
a<- working_dataset$gene_tissue_2
working_dataset$gene_tissue_2=NULL
working_dataset = as.data.frame(t(working_dataset))
colnames(working_dataset)<-a

origin_tissue = 'Adipose - Subcutaneous'
gene_tissue2 = paste0(gene_1, '_', origin_tissue)

tissue_set = unique(gsub(".*_","",tissue1$gene_tissue))
new_cols = met.brewer('Moreau', length(tissue_set))
names(new_cols) = unique(tissue_set)

origin = tissue1[,grepl(gene_tissue2, colnames(tissue1))]
target = tissue1[,!grepl(gene_tissue2, colnames(tissue1))]

full_cors = bicorAndPvalue(origin, target, use = 'p')
cor_table = reshape2::melt(full_cors$bicor)
cor_table$Var1=NULL
colnames(cor_table) = c('gene_tissue', 'bicor')
new_p = reshape2::melt(full_cors$p)

cor_table$pvalue = new_p$value[match(cor_table$gene_tissue, new_p$Var2)]
cor_table = na.omit(cor_table)
qest = qvalue(cor_table$pvalue, pi0 = 1)
cor_table$qvalue = qest$qvalues
cor_table$gene_symbol = gsub("\\_.*","",cor_table$gene_tissue)
cor_table$tissue = gsub(".*_","",cor_table$gene_tissue)
cor_table = cor_table[!is.na(cor_table$tissue),]
cor_table = cor_table[order(cor_table$pvalue, decreasing = F),]

new_table = annot[annot$pvalue<1e-4,]
new_table = new_table[!is.na(new_table$tissue_2),]
sig_set = new_table[!new_table$tissue_2==origin_tissue,]
row_length = ifelse(length(row.names(sig_set))>300, as.numeric(300), as.numeric(paste0(length(row.names(sig_set)))))
sig_set = sig_set[1:row_length,]
sig_set$gene_tissue_2 = paste0(sig_set$gene_symbol_2, '_', sig_set$tissue_2)

all_tog2 = working_dataset[, colnames(working_dataset) %in% sig_set$gene_tissue_2 | colnames(working_dataset)== gene_tissue2]
all_tog2[is.na(all_tog2)] = 0

colkey1 = as.data.frame(gsub(".*_","",colnames(all_tog2)))
colnames(colkey1) = 'tissue'

colkey1$cols = new_cols[match(colkey1$tissue, names(new_cols))]
colkey1<-na.omit(colkey1)
colkey1$cols<-as.numeric(colkey1$cols)
bics_map = bicorAndPvalue(all_tog2, all_tog2, use = 'p')

map1 = as.data.frame(bics_map$bicor)
map1 = reshape2::melt(as.matrix(map1))
map1$value[map1$value > 0.999999] <- 0
map2 = dcast(map1, Var1 ~ Var2, value.var = 'value')
row.names(map2) = map2$Var1
map2$Var1 = NULL

pdf(file = paste('Undirected network qvalLess 1e-3 - ', gene_1, '.pdf'))
qgraph(map2, minimum = 0.2, cut = 0.6, vsize = 3, color=colkey1$cols, legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=2, directed=F, labels = colnames(map2)) + ggtitle('')
dev.off()


net1 = mmhc(all_tog2) 
pdf(file = paste('Directed network qvalLess 1e-3 - ', gene_1, '.pdf'))
qgraph(net1, vsize = 3, legend = F, color=colkey1$cols, borders = TRUE, layout='spring', label.cex=2, directed=T, labels = colnames(map2)) + ggtitle('')
dev.off()

pdf(file = paste('Tissue Legend', gene_1, '.pdf'))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =names(new_cols), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = new_cols)
mtext("Tissue", at=0.2, cex=2)
dev.off()

## enrichment 
candid_gene = paste0(gene_1, '_', origin_tissue)
new_working1 = working_dataset[,grepl(candid_gene, colnames(working_dataset))]
new_working2 = working_dataset[,!grepl(candid_gene, colnames(working_dataset))]
full_cors = bicorAndPvalue(new_working1, new_working2, use = 'p')
cor_table = reshape2::melt(full_cors$bicor)
head(cor_table)
cor_table$Var1=NULL
colnames(cor_table) = c('gene_tissue', 'bicor')
new_p = reshape2::melt(full_cors$p)

cor_table$pvalue = new_p$value[match(cor_table$gene_tissue, new_p$Var2)]
cor_table = na.omit(cor_table)
qest = qvalue(cor_table$pvalue, pi0 = 1)
cor_table$qvalue = qest$qvalues
cor_table$gene_symbol = gsub("\\_.*","",cor_table$gene_tissue)
cor_table$tissue = gsub(".*_","",cor_table$gene_tissue)
cor_table = cor_table[!is.na(cor_table$tissue),]
tissue_col = met.brewer('Tara', length(unique(cor_table$tissue)))
names(tissue_col) = unique(cor_table$tissue)
cor_table$tissue_col = tissue_col[match(cor_table$tissue, names(tissue_col))]
cor_table = cor_table[order(cor_table$pvalue, decreasing = F),]

res1 = cor_table[cor_table$pvalue<0.01,]
res1 = na.omit(res1)
write.csv(res1, file = paste0('significant crosstissue enrichments with ',  candid_gene, ' -origin included.csv'), row.names = F)

sig_table = annot[annot$qvalue<0.1,]
sig_table$qcat =ifelse(sig_table$qvalue<0.01, 'q<0.01', 'q<0.1')
sig_table$qcat =ifelse(sig_table$qvalue<0.0001, 'q<0.0001', paste0(sig_table$qcat))
table(sig_table$qcat)

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
  dplyr::group_by(tissue_2, qcat) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))

pdf(file = paste0('significant crosstissue enrichments with ',  candid_gene, '.pdf'))
ggplot(binned_sig_prots, aes(x = "", y = freq, fill =tissue_2)) + 
  geom_bar(stat = "identity", width = 1, position = position_fill()) + blank_theme + theme(plot.title=element_text(size=25, face="bold")) +
  theme(axis.text.x=element_blank())+ scale_fill_manual(values=met.brewer('Moreau', length(unique(sig_table$tissue_2)))) +
  coord_polar(theta = "y") + 
  facet_wrap( ~ qcat)
dev.off()

tissue_list = binned_sig_prots[binned_sig_prots$qcat=='q<0.01',]
tissue_list = tissue_list[order(tissue_list$n, decreasing = T),]

select_tissue = tissue_list$tissue_1[1]
pp1 = annot[annot$qvalue<0.05,]
pp1 = pp1[pp1$tissue_2 %in% select_tissue,]
pp1_length = ifelse(length(row.names(pp1)) > 200, as.numeric(200), as.numeric(length(row.names(pp1))))
pp2 = pp1[1:pp1_length,]
gg1 = pp2$gene_symbol_2

setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs1 <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "Reactome_2022", "DSigDB")

enriched <- enrichr(gg1, dbs1)
names(enriched[1])
pdf(file = paste0(select_tissue, ' gene correlations with ', candid_gene, ' ', names(enriched[1]), '.pdf'))
plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()
pdf(file = paste0(select_tissue, ' gene correlations with ', candid_gene, ' ', names(enriched[2]), '.pdf'))
plotEnrich(enriched[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()
