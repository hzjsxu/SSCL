##TODO:
## analyze the paternal genes expression from hybrid embryo.

setwd('hybrid_RNA-seq')

library(rtracklayer)
library(GenomicFeatures)
library(GeneOverlap)
library(VennDiagram)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggplot2)
library(ggpubr)

get_time_point <- function(x, n_stages=1, expression_threshold=0.5){
  for (i in 1:(length(x) - n_stages)){
    time_point <- NA
    if (all(x[i:(i+n_stages)] > expression_threshold)){
      time_point <- as.numeric(gsub('S', '', names(x)[i]))
      break
    }
    else{
      next
    }
  }
  return(time_point)
}

## reading gtf file.
gtf_file <- "merge.XENTR_10.0_XENLA_9.2_GCF.gtf"
x_merge_gtf <- rtracklayer::import(gtf_file)
x_merge_gtf <- as.data.frame(x_merge_gtf)
x_merge_gene <- unique(x_merge_gtf[, c("gene_id", "gene_name")])

## reading gene expression file.
SAMPLES <- c('S3-1', 'S3-2',
             'S5-1', 'S5-2',
             'S7-1', 'S7-2',
             'S8-1', 'S8-2',
             'S9-1', 'S9-2',
             'S10-1', 'S10-2',
             'S11-1', 'S11-2',
             'S12-1', 'S12-2',
             'S13-1', 'S13-2',
             'S16-1', 'S16-2',
             'S24-1', 'S24-2',
             'S34-1', 'S34-2',
             'S41-1', 'S41-2')
XLXT <- x_merge_gene
gene_exp_files <- list.files(path = '.', '*.genes.results')
for (sample in SAMPLES){
  gene_exp_file <- paste0('XLXT-', sample, '.genes.results')
  cat(gene_exp_file, '\n')
  tmp <- read.table(gene_exp_file, sep = '\t', header = T)[, c('gene_id', 'TPM')]
  names(tmp) <- c('gene_id', sample)
  XLXT <- merge(XLXT, tmp, by='gene_id')
}

## get mean expression of two replicates.
df <- x_merge_gene
for (i in seq(1, length(SAMPLES), 2)){
  tmp <- data.frame(gene_id=XLXT$gene_id,
                    sample_name=(XLXT[, i+2] + XLXT[, i+3])/2)
  names(tmp)[2] <- strsplit(SAMPLES[i], '-')[[1]][1]
  df <- merge(df, tmp, by='gene_id')
}
XLXT <- df

## get gene expression derived from paternal (tropicalis) and maternal (laevis), respectively.
XLXT_tropicalis <- XLXT[1:28772, ]
XLXT_laevis <- XLXT[28773:62288, ]

## anchor genes (78).
anchor_genes <- read.table('Anchor_Genes.txt', header = F, sep = '\t', 
                           col.names = c('anchor_chr', 'anchor_start', 'anchor_end',
                                         'gene_chr', 'gene_start', 'gene_end', 'gene_name', 'gene_len', 'gene_strand'))
anchor_gene_exp <- XLXT_tropicalis[XLXT_tropicalis$gene_name %in% anchor_genes$gene_name, ]

## x.tropicalis protein-coding genes.
tropicalis_protein_coding_genes <- read.table('xv10.protein_coding.gene.bed', header = F, sep = '\t',
                                              col.names = c('gene_chr', 'gene_type', 'gene_start', 'gene_name', 'gene_len', 'gene_strand'))
tropicalis_protein_coding_genes_exp <- XLXT_tropicalis[XLXT_tropicalis$gene_name %in% tropicalis_protein_coding_genes$gene_name, ]

##### define function to get row-scaled expression matrix quickly and perform Z-score normalization.
rowScaleExpr <- function(expr_set){
  zscore_expr_set <- apply(expr_set, 1, scale)
  zscore_expr_set[is.na(zscore_expr_set)] <- 0
  zscore_expr_set <- as.data.frame(t(zscore_expr_set))
  names(zscore_expr_set) <- names(expr_set)
  zscore_expr_set
}

XLXT_exp <- data.frame(XLXT[, 3:ncol(XLXT)], row.names = XLXT$gene_id)
XLXT_exp_zscore <- rowScaleExpr(XLXT_exp)

XLXT_tropicalis_exp <- data.frame(XLXT_tropicalis[, 3:ncol(XLXT_tropicalis)], row.names = XLXT_tropicalis$gene_id)
XLXT_tropicalis_exp_zscore <- rowScaleExpr(XLXT_tropicalis_exp)

XLXT_laevis_exp <- data.frame(XLXT_laevis[, 3:ncol(XLXT_laevis)], row.names = XLXT_laevis$gene_id)
XLXT_laevis_exp_zscore <- rowScaleExpr(XLXT_laevis_exp)

anchor_gene_exp <- data.frame(anchor_gene_exp[, 3:ncol(anchor_gene_exp)], row.names = anchor_gene_exp$gene_id)
rownames(anchor_gene_exp) <- gsub('gene-', '', rownames(anchor_gene_exp))
anchor_gene_exp_zscore <- rowScaleExpr(anchor_gene_exp)

tropicalis_protein_coding_genes_exp <- data.frame(tropicalis_protein_coding_genes_exp[, 3:ncol(tropicalis_protein_coding_genes_exp)],
                                                  row.names = tropicalis_protein_coding_genes_exp$gene_id)
tropicalis_protein_coding_genes_exp_zscore <- rowScaleExpr(tropicalis_protein_coding_genes_exp)

#########  1. For all x.tropicalis genes.
### 1.1 log2 normalization.
df <- tropicalis_protein_coding_genes_exp[apply(tropicalis_protein_coding_genes_exp, 1, function(x){sd(x) != 0}),]
rownames(df) <- gsub('gene-', '', rownames(df))
pdf('all_protein_coding_genes.log2TPM.cluster_4.heatmap.pdf', height = 10, width=3)
tropicalis_protein_coding_genes_heatmap <- Heatmap(log(df + 1, 2),
                                                   border = NA,
                                                   name = 'log2(TPM+1)',
                                                   km = 4, cluster_rows = T,
                                                   cluster_columns = F,
                                                   show_row_names = F, row_names_gp = gpar(fontsize = 7),
                                                   show_column_names = T, column_names_gp = gpar(fontsize = 8),
                                                   clustering_distance_rows = 'pearson', clustering_method_rows = 'complete',
                                                   row_dend_reorder = T,
                                                   col = colorRamp2(seq(0, 5, length.out=100),
                                                                    colorRampPalette(c('#015f74', '#123634', '#e1e603'))(100)),
)
#sapply(row_order(tropicalis_protein_coding_genes_heatmap), length)
draw(tropicalis_protein_coding_genes_heatmap, column_title='19,522 expressed genes', column_title_gp=gpar(fontsize=14))
dev.off()

## 1.2 determine the best cluster: 3/4
mydata <- log(df + 1, 2)
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15){
  wss[i] <- sum(kmeans(mydata,centers=i)$withinss)
}
###这里的wss(within-cluster sum of squares)是组内平方和
df_variance <- data.frame(Clusters=1:15, Variance=wss)
pdf('all_genes.kmeans_cluster.pdf', height = 3, width = 3.5)
ggplot(data = df_variance, aes(x=Clusters, y=Variance))+
  geom_line()+
  geom_point()+
  scale_x_continuous(breaks = seq(1,15))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=12, face='bold', color='black'),
        plot.title = element_text(size=14, face='bold', color='black', hjust = 0.5),
  )
dev.off()

## 1.3 Heatmap according to reordered cluster.
cluster_info <- row_order(tropicalis_protein_coding_genes_heatmap)
df_new <- df[c(cluster_info$`1`, cluster_info$`2`, cluster_info$`3`, cluster_info$`4`), ]
df_unexpressed <- tropicalis_protein_coding_genes_exp[apply(tropicalis_protein_coding_genes_exp, 1, function(x){sd(x) == 0}),]
rownames(df_unexpressed) <- gsub('gene-', '', rownames(df_unexpressed))
df_new <- rbind(df_unexpressed, df_new)
rownames(df_new) <- gsub('gene-', '', rownames(df_new))

## reorder
df_new <- df[c(cluster_info$`4`,
               rev(cluster_info$`3`),
               cluster_info$`2`,
               cluster_info$`1`), ]
df_new <- rbind(df_new, df_unexpressed)
rownames(df_new) <- gsub('gene-', '', rownames(df_new))
pdf('protein_coding_genes.log2TPM.cluster_3.heatmap.pdf', height = 10, width=3)
tropicalis_protein_coding_genes_heatmap <- Heatmap(log(df_new + 1, 2),
                                            border = NA,
                                            row_split = rep(c('C1','C2','C3'), c(4572,5734,11443)),
                                            name = 'log2(TPM+1)',
                                            cluster_columns = F, cluster_rows = F,
                                            show_row_names = F, row_names_gp = gpar(fontsize = 7),
                                            show_column_names = T, column_names_gp = gpar(fontsize = 8),
                                            clustering_distance_rows = 'pearson', clustering_method_rows = 'complete',
                                            row_dend_reorder = T,
                                            col = colorRamp2(seq(0, 5, length.out=100),
                                                            colorRampPalette(c('#015f74', '#123634', '#e1e603'))(100)),
)
#sapply(row_order(tropicalis_protein_coding_genes_heatmap), length)
draw(tropicalis_protein_coding_genes_heatmap, column_title='21,749 protein-coding genes', column_title_gp=gpar(fontsize=14))
dev.off()


###### 2. For 78 x.tropicalis anchor genes.
### 2.1 zscore normalization.
pdf('anchor_genes.zscore_TPM.heatmap.pdf', height=8, width=4.5)
anchor_gene_heatmap <- Heatmap(anchor_gene_exp_zscore[apply(anchor_gene_exp_zscore, 1, function(x){sd(x) != 0}),],
                               border = TRUE,
                               name = 'Z-score TPM',
                               km = 3, cluster_rows = T,
                               cluster_columns = F,
                               show_row_names = T, row_names_gp = gpar(fontsize = 7),
                               show_column_names = T, column_names_gp = gpar(fontsize = 8),
                               clustering_distance_rows = 'spearman', clustering_method_rows = 'complete',
                               row_dend_reorder = T,
                               col = colorRamp2(seq(-3, 3, length.out=100),
                                                colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100))
)
draw(anchor_gene_heatmap, column_title='78 anchor genes', column_title_gp=gpar(fontsize=14))
dev.off()

### 2.2 log2 normalization.
pdf('anchor_genes.log2TPM.heatmap.pdf', height=8, width=4)
anchor_gene_heatmap <- Heatmap(log(anchor_gene_exp+1, 2),
                               border = TRUE,
                               name = 'log2(TPM+1)',
                               km = 3, cluster_rows = T,
                               cluster_columns = F,
                               show_row_names = T, row_names_gp = gpar(fontsize = 7),
                               show_column_names = T, column_names_gp = gpar(fontsize = 8),
                               clustering_distance_rows = 'euclidean', clustering_method_rows = 'complete',
                               row_dend_reorder = T,
                               col = colorRamp2(seq(0, 5, length.out=100),
                                                colorRampPalette(c('#015f74', '#123634', '#e1e603'))(100)))
draw(anchor_gene_heatmap, column_title='78 anchor genes', column_title_gp=gpar(fontsize=14))
dev.off()

### 2.3 anchor gene expression heatmap.
unexpressed_anchor_genes <- anchor_gene_exp[apply(anchor_gene_exp, 1, function(x){sd(x) == 0}),]
df_anchor_gene_exp_C1 <- anchor_gene_exp[rownames(anchor_gene_exp) %in% rownames(df_new[1:4572,]), ]
df_anchor_gene_exp_C2 <- anchor_gene_exp[rownames(anchor_gene_exp) %in% rownames(df_new[(4572+1):(4572+5734),]), ]
df_anchor_gene_exp_C3 <- anchor_gene_exp[rownames(anchor_gene_exp) %in% rownames(df_new[(4572+5734+1):(4572+5734+11443),]), ]
df_anchor_gene_exp <- rbind(df_anchor_gene_exp_C1) %>%
  rbind(df_anchor_gene_exp_C2) %>%
  rbind(df_anchor_gene_exp_C3) %>%
  rbind(unexpressed_anchor_genes)

pdf('anchor_genes.log2TPM.cluster_3.heatmap.pdf', height=8, width=4)
anchor_gene_heatmap <- Heatmap(log(df_anchor_gene_exp+1, 2),
                               border = TRUE,
                               row_split = rep(c('C1','C2','C3'), c(10,19,49)),
                               name = 'log2(TPM+1)',
                               cluster_columns = F, cluster_rows = F,
                               show_row_names = T, row_names_gp = gpar(fontsize = 7),
                               show_column_names = T, column_names_gp = gpar(fontsize = 8),
                               clustering_distance_rows = 'euclidean', clustering_method_rows = 'complete',
                               row_dend_reorder = T,
                               col = colorRamp2(seq(0, 5, length.out=100),
                                                colorRampPalette(c('#015f74', '#123634', '#e1e603'))(100)))
draw(anchor_gene_heatmap, column_title='78 anchor genes', column_title_gp=gpar(fontsize=14))
dev.off()



#### 3. percentage of three clusters.
df1 <- data.frame(cluster=rep(c('C1','C2','C3'), 2),
                  gene_type=rep(c('all', 'anchor_gene'), each=3),
                  percentage=c(21.02, 26.37, 52.61,
                               12.82, 24.36, 62.82))
df1$gene_type <- factor(df1$gene_type, levels = c('all', 'anchor_gene'))

pdf('./Three_Cluster_Genes.Percentage.barplot.pdf', width=6, height=4)
ggplot(data = df1, aes(x = cluster, y = percentage, fill=gene_type))+
  geom_bar(position=position_dodge(0.6),stat="identity",width = 0.6)+
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, 20))+
  scale_fill_manual(values = c('#bdbdbd', '#6a7de4'))+
  geom_text(aes(x = cluster, y = percentage, label=percentage), size=4, position=position_dodge(0.6),stat="identity",vjust=-0.4)+
  labs(x='', y='Percentage (%)')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face='plain'),
        axis.text.x = element_text(hjust = 0.5, vjust=0.5, size = 8, face='plain', color='black'),
        axis.text.y = element_text(hjust = 0.5, size = 8, face='plain', color='black'),
        axis.title = element_text(size = 10, color='black', face='plain'))
dev.off()

#### permutation test
all_Cluster_Genes <-data.frame(genename=rownames(df_new),
                               cluster=rep(c('C1','C2','C3'), c(4572, 5734, 11443)))
anchor_Cluster_Genes <- data.frame(genename=rownames(df_anchor_gene_exp),
                                   cluster=rep(c('C1','C2','C3'), c(10,19,49)))
Superloop_Cluster_Genes <- data.frame(genename=rownames(df_loop_cluster_gene_exp),
                                      cluster=rep(c('C1','C2','C3','C4'), c(223, 751, 1334, 2000)))
df_permutation <- data.frame(Times=character(),
                             C1=numeric(),
                             C2=numeric(),
                             C3=numeric())
n_anchor_C1_genes <- nrow(anchor_Cluster_Genes[anchor_Cluster_Genes$cluster == 'C1', ])
n_anchor_C2_genes <- nrow(anchor_Cluster_Genes[anchor_Cluster_Genes$cluster == 'C2', ])
n_anchor_C3_genes <- nrow(anchor_Cluster_Genes[anchor_Cluster_Genes$cluster == 'C3', ])

for (i in seq(1,1000)){
  time=paste0('Permutation_', i)
  cat('Permutation ', i, '...\n')
  random_C1_genes <- sample(all_Cluster_Genes$genename, n_anchor_C1_genes)
  random_C2_genes <- sample(all_Cluster_Genes$genename, n_anchor_C2_genes)
  random_C3_genes <- sample(all_Cluster_Genes$genename, n_anchor_C3_genes)
  n_random_C1_anchor_genes <- round(length(random_C1_genes[random_C1_genes %in% all_Cluster_Genes[all_Cluster_Genes$cluster == 'C1', ]$genename])/n_anchor_C1_genes,2)
  n_random_C2_anchor_genes <- round(length(random_C2_genes[random_C2_genes %in% all_Cluster_Genes[all_Cluster_Genes$cluster == 'C2', ]$genename])/n_anchor_C2_genes,2)
  n_random_C3_anchor_genes <- round(length(random_C3_genes[random_C3_genes %in% all_Cluster_Genes[all_Cluster_Genes$cluster == 'C3', ]$genename])/n_anchor_C3_genes,2)
  df_permutation[i, ] <- c(time, n_random_C1_anchor_genes, n_random_C2_anchor_genes, n_random_C3_anchor_genes)
}
df_permutation$C1 <- as.numeric(df_permutation$C1)
df_permutation$C2 <- as.numeric(df_permutation$C2)
df_permutation$C3 <- as.numeric(df_permutation$C3)
# df_permutation$C4 <- as.numeric(df_permutation$C4)

t.test(df_permutation$C1, mu = 0.1282, alternative = 'two.sided') # P<2.2×10-16
t.test(df_permutation$C2, mu = 0.2436, alternative = 'two.sided') # P=3.622e-07
t.test(df_permutation$C3, mu = 0.6282, alternative = 'two.sided') # P<2.2×10-16

## 4. stat expressed genes in each stage.
########## 4.1 For X.tropicalis genes.
df <- tropicalis_protein_coding_genes_exp[apply(tropicalis_protein_coding_genes_exp, 1, function(x){sd(x) != 0}),]
df_num_XT_expressed_genes <- as.data.frame(apply(df, 2, function(x) length(x[x>0.5])))  ## TPM > 0.5
df_num_XT_expressed_genes <- data.frame(Stage=rownames(df_num_XT_expressed_genes),
                                        Number=df_num_XT_expressed_genes$`apply(df, 2, function(x) length(x[x > 0.5]))`)
df_num_XT_expressed_genes$Stage <- factor(df_num_XT_expressed_genes$Stage, levels = rev(df_num_XT_expressed_genes$Stage))
pdf('Numbers.expressed.genes.in.each.stage.pdf', width = 6, height = 4)
ggplot(df_num_XT_expressed_genes, aes(x = Stage, y = Number))+
  geom_bar(stat = 'identity', width = 0.6, fill='#4574b8', color='black')+
  scale_y_continuous(limits = c(0, 20000), breaks = seq(0, 20000, 5000))+
  geom_text(aes(x = Stage, y = Number, label=Number), size=4, hjust=-0.1)+
  coord_flip()+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face='bold'),
        axis.text.x = element_text(hjust = 0.5, vjust=0.5, size = 8, face='bold', color='black'),
        axis.text.y = element_text(hjust = 0.5, size = 8, face='bold', color='black'),
        axis.title = element_text(size = 10, color='black', face='bold'),
        legend.position = 'none')
dev.off()

########## 4.2 For X.laevis genes.
df_laevis <- XLXT_laevis_exp[apply(XLXT_laevis_exp, 1, function(x){sd(x) != 0}),]
df_num_XL_expressed_genes <- as.data.frame(apply(df_laevis, 2, function(x) length(x[x>0.5])))  ## TPM > 0.5
df_num_XL_expressed_genes <- data.frame(Stage=rownames(df_num_XL_expressed_genes),
                                        Number=df_num_XL_expressed_genes$`apply(df_laevis, 2, function(x) length(x[x > 0.5]))`)
df_num_XL_expressed_genes$Stage <- factor(df_num_XL_expressed_genes$Stage, levels = rev(df_num_XL_expressed_genes$Stage))
pdf('Numbers.XL.expressed.genes.in.each.stage.pdf', width = 6, height = 4)
ggplot(df_num_XL_expressed_genes, aes(x = Stage, y = Number))+
  geom_bar(stat = 'identity', width = 0.6, fill='#de773f', color='black')+
  scale_y_continuous(limits = c(0, 25000), breaks = seq(0, 25000, 5000))+
  geom_text(aes(x = Stage, y = Number, label=Number), size=4, hjust=-0.1)+
  coord_flip()+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face='bold'),
        axis.text.x = element_text(hjust = 0.5, vjust=0.5, size = 8, face='bold', color='black'),
        axis.text.y = element_text(hjust = 0.5, size = 8, face='bold', color='black'),
        axis.title = element_text(size = 10, color='black', face='bold'),
        legend.position = 'none')
dev.off()

######### 4.3 gene expression for x.tropicalis and x.laevis.
df_num_paternal_expressed_genes <- data.frame(Stage=df_num_XT_expressed_genes$Stage,
                                              Number=c(df_num_XT_expressed_genes$Number,
                                                       df_num_XL_expressed_genes$Number),
                                              Source=c(rep('XT', 13),
                                                       rep('XL', 13)))
df_num_paternal_expressed_genes$Stage <- factor(df_num_paternal_expressed_genes$Stage,
                                                levels = df_num_XT_expressed_genes$Stage)
pdf('./Numbers.paternal.expressed.genes.in.each.stage.pdf', width = 10, height = 5)
ggplot(df_num_paternal_expressed_genes, aes(x = Stage, y = Number, fill=Source))+
  geom_bar(position = 'dodge', stat = 'identity', width = 0.8)+
  scale_y_continuous(limits = c(0, 25000), breaks = seq(0, 25000, 5000))+
  geom_text(aes(x = Stage, y = Number, label=Number), size=4, vjust=-0.2, hjust=0.5, position=position_dodge(width = 0.8))+
  scale_fill_manual(values=c('#4574b8', '#de773f'))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face='bold'),
        axis.text.x = element_text(hjust = 0.5, vjust=0.5, size = 8, face='bold', color='black'),
        axis.text.y = element_text(hjust = 0.5, size = 8, face='bold', color='black'),
        axis.title = element_text(size = 10, color='black', face='bold'),
  )
dev.off()

### 5. C1,C2,C3 expression time boxplot.
all_protein_coding_gene_exp_stage <- data.frame(hpf=apply(df_new, 1, get_time_point))
anchor_gene_exp_stage <- data.frame(hpf=apply(df_anchor_gene_exp, 1, get_time_point))
time_data <- data.frame(class=c(rep('C1', 4572),
                                rep('C2', 5734),
                                rep('C3', 11443)),
                        stage=all_protein_coding_gene_exp_stage$hpf)
time_data$class <- factor(time_data$class, levels = c('C3', 'C2', 'C1'))
comparions <- list(c('C1', 'C2'),
                   c('C2', 'C3'),
                   c('C1', 'C3'))
pdf('exp_stage.boxplot.pdf', width = 8, height = 3)
ggplot(data = time_data, aes(x=class, y=stage, fill=class))+
  stat_boxplot(geom='errorbar', width=0.3,)+
  geom_boxplot(width=0.5, outlier.colour = NA)+
  stat_compare_means(comparisons = comparions) +
  scale_fill_manual(values=c('#56719f', '#af5a5d', '#c68560', '#009ad6'))+
  labs(title = "expression time of C1-C3 genes", x='', y='Stage')+
  coord_flip()+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face='bold'),
        axis.text.x = element_text(hjust = 0.5, vjust=0.5, size = 10, face='plain', color='black'),
        axis.text.y = element_text(hjust = 0.5, size = 10, face='plain', color='black'),
        axis.title = element_text(size = 12, color='black', face='plain'),
        legend.position = 'none')
dev.off()