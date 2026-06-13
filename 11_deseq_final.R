library(ggplot2)
library(DESeq2)
library(umap)
library(ggpubr)
library(viridis)
library(ggrepel)
library(ggVennDiagram)
library(pheatmap)

################ prepare input data

set.seed(12345)

# load raw counts for protein-coding genes
gene_counts <- read.delim(file = "9_count_reads/gene_counts_pc.tsv",
                          header = T, sep = "\t",
                          stringsAsFactors = F,
                          check.names = F)
head(gene_counts)

# separate metadata and counts
counts_tab <- gene_counts[, grepl(pattern = "_lt2", x = colnames(gene_counts))]
rownames(counts_tab) <- paste(gene_counts$locus_tag, gene_counts$gene_name, sep = "_")
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


# note: deseq needs column names here to match rownames of count table.
# they must also be in the same order
rownames(cond) <- cond$sample
table(colnames(counts_tab)==rownames(cond))
table(rownames(cond)==colnames(counts_tab))

################ get variance stabilized estimate

# remove samples with lots of zeros
colSums(counts_tab)
zeroprop <- apply(X = counts_tab, MARGIN = 2, FUN = function(x) length(which(x == 0))/length(x))
zeroprop
length(counts_tab$`1_trial1_day0_malic-acid_lt2`[which(counts_tab$`1_trial1_day0_malic-acid_lt2`==0)])/nrow(counts_tab)
# remove samples with >= 0.90
bad <- zeroprop[which(zeroprop >= 0.90)]
bad

counts_tab <- counts_tab[,which(!(colnames(counts_tab)%in%names(bad)))]
head(counts_tab)

cond <- cond[which(cond$sample%in%colnames(counts_tab)),]

################ differential expression testing

# create DESeq object from count matrix, this time with group in our design
dds <- DESeqDataSetFromMatrix(countData = counts_tab,
                              colData = cond,
                              design= ~ group) 

# run DESeq
dds <- DESeq(dds)

################ malic-acid vs control

# test malic-acid vs control
res <- results(object = dds, c("group", "malic-acid", "control"))

# reorder based on p-value
res <- res[order(res$pvalue),] 

# add name column with gene IDs
res$name <- rownames(res)
head(res)

# add siggene column
# adjusted p-value < 0.05 & log2fc >= 1 == Upregulated
# adjusted p-value < 0.05 & log2fc <= -1 == Downregulated
# Everything else == Not significant
res$siggene <- ifelse(test = (res$padj<0.05 & res$log2FoldChange>=1.0), yes = "Upregulated",
                    no = ifelse(test = (res$padj<0.05 & res$log2FoldChange<=-1.0), yes = "Downregulated",
                                no = "Not significant"))

table(res$siggene)

# get quantiles of -log10(p-values)
quantile(x = -log10(res$pvalue),
         probs = c(0.75, 0.90, 0.95, 0.99), na.rm = T)

# plot histogram of -log10(p-values)
ggplot(data = as.data.frame(res), mapping = aes(x = -log10(pvalue))) + 
  geom_histogram() +
  theme_bw()

# get most significant p-values
supersig <- rownames(res)[which(-log10(res$pvalue)>=5.34)]
supersig

# plot volcano plot
pdf(file = "11_deseq_final/volcano_lt2.pdf", width = 8, height = 8)
ggplot(data = as.data.frame(res), aes(x = log2FoldChange, y = -log10(pvalue), color = siggene)) + 
  geom_point() + 
  theme_bw() +
  scale_color_manual(values = c("blue3", "gray54", "red3")) +
  geom_text_repel(aes(label = ifelse(name%in%supersig, name, "")))
dev.off()

# save results to TSV
write.table(x = res, file = "11_deseq_final/deseq_lt2.tsv",
            append = F, quote = F,
            sep = "\t", row.names = T,
            col.names = T)

################ heatmap of most significant genes

# get the 50 most significant genes
supertop <- res
supertop <- supertop[order(supertop$pvalue, decreasing = F),]
head(supertop)

all_topgenes <- supertop[1:50,]
dim(all_topgenes)
head(all_topgenes)
summary(all_topgenes$padj)

# scale and center raw counts
counts.scaled <- scale(x = t(counts_tab), center = T, scale = T)
head(counts.scaled[,1:5])

# plot heatmap of most significant genes
# note that batch effects not removed; interpret with caution
row.metadata <- unlist(lapply(strsplit(x = rownames(counts.scaled), split = "_"), "[[", 4))
row.metadata <- as.data.frame(gsub(pattern = "-", replacement = "_", x = row.metadata))
rownames(row.metadata) <- rownames(counts.scaled)
colnames(row.metadata) <- "Group"
row.metadata

ann_colors <- list(Group = c(control="#66CCEE", malic_acid="#AA3377"))
ann_colors


# save heatmap as PDF
pdf(file = "11_deseq_final/heatmap_top50genes.pdf", width = 11, height = 8)
pheatmap(mat = counts.scaled[,which(colnames(counts.scaled)%in%all_topgenes$name)],
         fontsize_col = 5, fontsize_row = 8,
         annotation_row = row.metadata,
         annotation_colors = ann_colors)
dev.off()
