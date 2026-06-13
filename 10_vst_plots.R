library(ggplot2)
library(ggpubr)
library(DESeq2)
library(umap)
library(dendextend)
library(pvclust)
library(corrplot)
library(ggrepel)
library(viridis)
library(khroma)
library(gridExtra)

# set seed
set.seed(12345)

################ prepare input data

# load raw counts for protein-coding genes
gene_counts <- read.delim(file = "9_count_reads/gene_counts_pc.tsv",
                        header = T, sep = "\t",
                        stringsAsFactors = F,
                        check.names = F)
head(gene_counts)

# separate metadata and counts
counts_tab <- gene_counts[, grepl(pattern = "_lt2", x = colnames(gene_counts))]
rownames(counts_tab) <- gene_counts$feature_id
head(counts_tab)

meta_tab <- gene_counts[, !(grepl(pattern = "_lt2", x = colnames(gene_counts)))]
head(meta_tab)

# make data frame of conditions
cond <- data.frame(sample = colnames(counts_tab))
cond

# make trial column
cond
cond$trial <- unlist(lapply(strsplit(x = cond$sample, split = "_"), "[[", 2))
cond
table(cond$trial)

# make group column
cond$group <- unlist(lapply(strsplit(x = cond$sample, split = "_"), "[[", 4))
cond
table(cond$group)

# make day column
cond$day <- unlist(lapply(strsplit(x = cond$sample, split = "_"), "[[", 3))
cond
table(cond$day)

# make sequencing depth column
table(colnames(counts_tab)%in%cond$sample)
table(cond$sample%in%colnames(counts_tab))
table(colnames(counts_tab)==cond$sample)
table(cond$sample==colnames(counts_tab))

cond$seqdepth <- colSums(counts_tab)
cond


# note: deseq need column names here to match rownames of count table.
# they must also be in the same order
rownames(cond) <- cond$sample
table(colnames(counts_tab)==rownames(cond))
table(rownames(cond)==colnames(counts_tab))

################ get variance stabilized estimate

# remove those with lots of zeros
colSums(counts_tab)
zeroprop <- apply(X = counts_tab, MARGIN = 2, FUN = function(x) length(which(x == 0))/length(x))
zeroprop
length(counts_tab$`1_trial1_day0_malic-acid_lt2`[which(counts_tab$`1_trial1_day0_malic-acid_lt2`==0)])/nrow(counts_tab)
# remove those with >= 0.90
bad <- zeroprop[which(zeroprop >= 0.90)]
bad

counts_tab <- counts_tab[,which(!(colnames(counts_tab)%in%names(bad)))]
head(counts_tab)

cond <- cond[which(cond$sample%in%colnames(counts_tab)),]

# create DESeq data set object
dds <- DESeqDataSetFromMatrix(countData = counts_tab,
                              colData = cond,
                              design= ~ 1) # initialize with no design

# run DESeq 
dds <- DESeq(dds)

# apply variance stabilizing transformation (vst)
vsd <- varianceStabilizingTransformation(dds)

# store vst-transformed counts as variable foo
foo <- vsd@assays@data[[1]]

# sanity check: column names of vst-transformed counts match sample order in conditions data frame
table(colnames(foo)==cond$sample)
table(cond$sample==colnames(foo))

################ hierarchical clustering of vst-transformed counts

# hierarchical clustering of vst-transformed counts
h <- pvclust(data = foo, method.hclust = "average", method.dist = "euclidean", nboot = 1000)
# plot basic dendrogram
plot(h)

# make it nicer
pdf(file = "10_vst_plots_expanded/sample_hclust.pdf", width = 8.5, height = 11)
par(mar = c(15, 3, 1, 1))
h.dend <-
  as.dendrogram(h) %>%
  set("branches_lwd", 3) %>% 
  pvclust_show_signif(h) %>% 
  plot(horiz = F)
h %>% text()
dev.off()

################ PCA of of vst-transformed counts
pca <- prcomp(t(foo), scale = F, center = T)
pca3col <- as.data.frame(cbind(pca$x[,1],pca$x[,2],pca$x[,3]))

table(rownames(pca3col)==cond$sample)
table(cond$sample==rownames(pca3col))
cond$PC1 <- pca3col$V1
cond$PC2 <- pca3col$V2
cond$PC3 <- pca3col$V3

# plot PCA, color by group
pca.group.hull <- ggplot(data = cond, aes(x = PC1, y = PC2, color = group)) + 
  geom_point(size = 3) +
  stat_chull(aes(color = group, fill = group), geom = "polygon", alpha = 0.1) + 
  theme_bw() +
  scale_fill_manual(values = c(control="#66CCEE", `malic-acid`="#AA3377")) + 
  scale_color_manual(values = c(control="#66CCEE", `malic-acid`="#AA3377"))
pca.group.hull

pca.group <- ggplot(data = cond, aes(x = PC1, y = PC2, size = PC3, color = group)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c(control="#66CCEE", `malic-acid`="#AA3377"))
pca.group

pca.group.text <- ggplot(data = cond, aes(x = PC1, y = PC2, size = PC3, color = group, label = sample)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c(control="#66CCEE", `malic-acid`="#AA3377")) +
  geom_text_repel()
pca.group.text

pdf(file = "10_vst_plots_expanded/pca_group.pdf", width = 11, height = 5)
grid.arrange(pca.group.hull, pca.group, pca.group.text, nrow = 1)
dev.off()

# plot PCA, color by trial
pca.trial.hull <- ggplot(data = cond, aes(x = PC1, y = PC2, color = trial)) + 
  geom_point(size = 3) +
  stat_chull(aes(color = trial, fill = trial), geom = "polygon", alpha = 0.1) + 
  theme_bw() 
  #scale_fill_manual(values = light(8)) + 
  #scale_color_manual(values = light(8))
pca.trial.hull

pca.trial <- ggplot(data = cond, aes(x = PC1, y = PC2, size = PC3, color = trial)) + 
  geom_point() +
  theme_bw() 
  #scale_color_manual(values = light(8))
pca.trial

pdf(file = "10_vst_plots_expanded/pca_trial.pdf", width = 11, height = 5)
grid.arrange(pca.trial.hull, pca.trial, nrow = 1)
dev.off()

# sequencing depth
pdf(file = "10_vst_plots_expanded/pca_seqdepth.pdf", width = 11, height = 5)
ggplot(data = cond, aes(x = PC1, y = PC2, size = PC3, color = log10(seqdepth), label = sample)) + 
  geom_point() +
  theme_bw() + 
  scale_color_viridis() + 
  geom_text_repel()
dev.off()
