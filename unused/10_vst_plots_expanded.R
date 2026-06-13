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
# par(mar = c(bottom, left, top, right))
pdf(file = "10_vst_plots_expanded/sample_hclust.pdf", width = 8.5, height = 11)
par(mar = c(15, 3, 1, 1))
h.dend <-
  as.dendrogram(h) %>%
  set("branches_lwd", 3) %>% 
  pvclust_show_signif(h) %>% 
  plot(horiz = F)
h %>% text()
dev.off()

# test if linkage method affects dendrogram
dend1 <- as.dendrogram(hclust(d = dist(t(foo)), method = "average"))
dend2 <- as.dendrogram(hclust(d = dist(t(foo)), method = "complete"))
dend3 <- as.dendrogram(hclust(d = dist(t(foo)), method = "centroid"))
dend4 <- as.dendrogram(hclust(d = dist(t(foo)), method = "median"))
dend5 <- as.dendrogram(hclust(d = dist(t(foo)), method = "single"))
dend6 <- as.dendrogram(hclust(d = dist(t(foo)), method = "ward.D"))
dend7 <- as.dendrogram(hclust(d = dist(t(foo)), method = "ward.D2"))
dend8 <- as.dendrogram(hclust(d = dist(t(foo)), method = "mcquitty"))

dend1to8 <- dendlist("Average" = dend1, "Complete" = dend2,
                     "Centroid" = dend3, "Median" = dend4,
                     "Single" = dend5, "Ward" = dend6,
                     "Ward.D2" = dend7, "McQuitty" = dend8)

pdf(file = "10_vst_plots_expanded/sample_hclust_corrplot.pdf", width = 8.5, height = 8.5)
par(mfrow=c(2,2))
corrplot(cor.dendlist(dend1to8), method = "pie", type = "lower")
corrplot(cor.dendlist(dend1to8), method = "ellipse", type = "lower")
corrplot(cor.dendlist(dend1to8), method = "number", type = "lower")
corrplot(cor.dendlist(dend1to8), method = "color", type = "lower")
dev.off()

# compare dend1 (average linkage) vs dend6 (Ward)
# Ward is the least similar to average
# calculate Baker’s Gamma Index (similarity measure between dendrograms)
# see https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html#correlation-measures
cor_bakers_gamma(dend1, dend6)

# perform permutation test to calculate statistical significance of index
# look at distribution of Baker’s Gamma Index under the null hypothesis (assuming fixed tree topologies)
# Here are the results when the compared tree is itself (after shuffling its own labels)...
# ...and when comparing tree 1 to the shuffled tree 2
nsamp <- 100
cor_bakers_gamma_results <- numeric(nsamp)
dend_mixed <- dend1
for(i in 1:nsamp) {
  dend_mixed <- sample.dendrogram(dend_mixed, replace = F)
  cor_bakers_gamma_results[i] <- cor_bakers_gamma(dend1, dend_mixed)
}

the_cor <- cor_bakers_gamma(dend1, dend1)
the_cor2 <- cor_bakers_gamma(dend1, dend6)

# plot results
pdf(file = "10_vst_plots_expanded/bakers_gamma_permutation_averageVSward.pdf", width = 8.5, height = 11)
plot(density(cor_bakers_gamma_results),
     main = "Baker's gamma distribution under H0",
     xlim = c(-1,1))
abline(v = 0, lty = 2)
abline(v = the_cor, lty = 2, col = 2)
abline(v = the_cor2, lty = 2, col = 4)
legend("topleft", legend = c("cor (self)", "cor2 (other)"), fill = c(2,4))
round(sum(the_cor2 < cor_bakers_gamma_results)/ nsamp, 4)
title(sub = paste("One sided p-value:",
                  "cor =",  round(sum(the_cor < cor_bakers_gamma_results)/ nsamp, 4),
                  " ; cor2 =",  round(sum(the_cor2 < cor_bakers_gamma_results)/ nsamp, 4)
))
dev.off()

# build bootstrap confidence interval, using sample.dendrogram, for the correlation
dend1_labels <- labels(dend1)
dend6_labels <- labels(dend6)
cor_bakers_gamma_results <- numeric(nsamp)
for(i in 1:nsamp) {
  sampled_labels <- sample(dend1_labels, replace = T)
  # members needs to be fixed since it will be later used in nleaves
  dend_mixed1 <- sample.dendrogram(dend1, 
                                   dend_labels=dend1_labels,
                                   fix_members=T, fix_order=T, fix_midpoint=F,
                                   replace = T, sampled_labels=sampled_labels
  )
  dend_mixed6 <- sample.dendrogram(dend6, dend_labels=dend6_labels,
                                   fix_members=T, fix_order=T, fix_midpoint=F,
                                   replace = T, sampled_labels=sampled_labels
  )                                    
  cor_bakers_gamma_results[i] <- cor_bakers_gamma(dend_mixed1, dend_mixed6, warn = F)
}

# plot the tanglegram
pdf(file = "10_vst_plots_expanded/tanglegram_averageVSward.pdf", width = 11, height = 8.5)
tanglegram(dend1, dend6)
dev.off()

# sanity check
cor_bakers_gamma(dend_mixed1, dend_mixed6, warn = F)

# get confidence interval
CI95 <- quantile(cor_bakers_gamma_results, probs=c(.025,.975))
CI95

# plot results 
# bootstrap sampling can do weird things with small trees
# interpret with caution!
pdf(file = "10_vst_plots_expanded/bakers_gamma_permutation_averageVSward_bootstrap.pdf", width = 8.5, height = 11)
par(mfrow = c(1,1))
plot(density(cor_bakers_gamma_results),
     main = "Baker's gamma bootstrap distribution",
     xlim = c(-1,1))
abline(v = CI95, lty = 2, col = 3)
abline(v = cor_bakers_gamma(dend1, dend6), lty = 2, col = 2)
legend("topleft", legend =c("95% CI", "Baker's Gamma Index"), fill = c(3,2))
dev.off()

# get cophenetic distance
cor_cophenetic(dend1, dend6)

# make Bk plot for Fowlkes-Mallows Index
pdf(file = "10_vst_plots_expanded/Bk_plot_averageVSward.pdf", width = 8.5, height = 8.5)
Bk_plot(dend1, dend6, main = "Bk plot \n(based on dendrograms)")
dev.off()

################ UMAP of of vst-transformed counts

# run UMAP on vst-transformed counts
out <- umap(t(foo), n_neighbors = 2)

# add UMAP results to condition data frame
table(rownames(out$layout)==cond$sample)
table(cond$sample==rownames(out$layout))
cond$umap_x <- out$layout[,1]
cond$umap_y <- out$layout[,2]

# plot UMAP results, colored by group
umap.group.hull <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = group)) + 
  geom_point(size = 3) +
  stat_chull(aes(color = group, fill = group), geom = "polygon", alpha = 0.1) + 
  theme_bw() +
  scale_fill_manual(values = c(control="#66CCEE", `malic-acid`="#AA3377")) + 
  scale_color_manual(values = c(control="#66CCEE", `malic-acid`="#AA3377"))
umap.group.hull

umap.group <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = group)) + 
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = c(control="#66CCEE", `malic-acid`="#AA3377"))
umap.group

umap.group.text <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = group, label = sample)) + 
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = c(control="#66CCEE", `malic-acid`="#AA3377")) +
  geom_text_repel()
umap.group.text

pdf(file = "10_vst_plots_expanded/umap_group.pdf", width = 11, height = 5)
grid.arrange(umap.group.hull, umap.group, umap.group.text, nrow = 1)
dev.off()

# plot UMAP results, colored by trial
light <- color("light")

cond$trial <- as.character(cond$trial)

umap.trial.hull <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = trial)) + 
  geom_point(size = 3) +
  stat_chull(aes(color = trial, fill = trial), geom = "polygon", alpha = 0.1) + 
  theme_bw() + 
  scale_fill_manual(values = light(8)) + 
  scale_color_manual(values = light(8))
umap.trial.hull

umap.trial <- ggplot(data = cond, aes(x = umap_x, y = umap_y, color = trial)) + 
  geom_point(size = 3) +
  theme_bw() +
  scale_fill_manual(values = light(8)) + 
  scale_color_manual(values = light(8))
umap.trial

pdf(file = "10_vst_plots_expanded/umap_trial.pdf", width = 11, height = 5)
grid.arrange(umap.trial.hull, umap.trial, nrow = 1)
dev.off()

# sequencing depth
pdf(file = "10_vst_plots_expanded/umap_seqdepth.pdf", width = 8.5, height = 5)
ggplot(data = cond, aes(x = umap_x, y = umap_y, color = log10(seqdepth), label = sample)) + 
  geom_point(size = 3) +
  theme_bw() + 
  scale_color_viridis() + 
  geom_text_repel()
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

