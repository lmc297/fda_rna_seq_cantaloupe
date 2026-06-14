library(ggplot2)
library(clusterProfiler)

# uses clusterProfiler: https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html

###################### set seed

set.seed(42)

###################### define functions: kegg pathway enrichment

# define run.kegg.pathway function, to perform kegg pathway enrichment
run.kegg.pathway <- function(genes.og, updown, fc, organism, p.cutoff, p.correction, treatment){
  
  # get significant genes, which are < user-supplied p-value cutoff
  genes <- genes.og[which(genes.og$padj < p.cutoff),]
  
  # if users want upregulated genes:
  # take significant genes with fold change >= user-supplied fold-change cutoff
  if (updown == "up"){
    sig.genes <- rownames(genes[which(genes$log2FoldChange >= fc),])
  }
  
  # if users want downregulated genes:
  # take significant genes with fold change <= user-supplied fold-change cutoff
  else if (updown == "down"){
    sig.genes <- rownames(genes[which(genes$log2FoldChange <= fc),])
  }
  
  # throw an error if user doesn't specify up or down
  else{
    stop("Please select 'up' or 'down' to query over-represented pathways among
              up- andn down-regulated genes, respectively.")
  }
  
  # KEGG pathway over-representation analysis
  # use enricher
  kk <- enrichKEGG(gene = sig.genes,
                   organism = organism,
                   pvalueCutoff = p.cutoff,
                   pAdjustMethod = p.correction,
                   minGSSize = 1, maxGSSize = 1e9)
  
  # KEGG pathway gene set enrichment analysis
  # rank genes by log2 fold change from highest to lowest
  genes.og <- genes.og[order(genes.og$log2FoldChange, decreasing = T),]
  allgenes.ranked <- genes.og$log2FoldChange
  names(allgenes.ranked) <- rownames(genes.og)
  print(allgenes.ranked)[1:5]
  
  # use GSEA
  kk2 <- gseKEGG(geneList = allgenes.ranked,
                 organism = organism,
                 minGSSize = 1,
                 maxGSSize = 1e9,
                 eps = 0,
                 pvalueCutoff = p.cutoff,
                 pAdjustMethod = p.correction,
                 verbose = F)
  
  # get results with adjusted p<0.05
  if (!(is.null(kk))){
  kk.final <- kk@result[which(kk@result$p.adjust < 0.05),]}
  else{
    kk.final <- data.frame()}
  if (!(is.null(kk2))){
  kk2.final <- kk2@result[which(kk2@result$p.adjust < 0.05),]}
  else{
    kk2.final <- data.frame()}
  
  # add treatment columns to final data frames
  kk.final$treatment <- rep(treatment, nrow(kk.final))
  kk2.final$treatment <- rep(treatment, nrow(kk2.final))
  
  # return results for both analyses
  return(list("pathway.overrep" = kk.final,
              "pathway.gsea" = kk2.final,
              "OGpathway.overrep" = kk@result,
              "OGpathway.gsea" = kk2@result))
}

###################### define functions: kegg module enrichment

# define run.kegg.module function to perform kegg module enrichment
run.kegg.module <- function(genes.og, updown, fc, organism, p.cutoff, p.correction, treatment){
  
  # get genes that have p-value < user-specified cutoff
  genes <- genes.og[which(genes.og$padj < p.cutoff),]
  
  # get significant up or down-regulated genes (or exit if up or down not specified)
  if (updown == "up"){
    sig.genes <- rownames(genes[which(genes$log2FoldChange >= fc),])
  }
  else if (updown == "down"){
    sig.genes <- rownames(genes[which(genes$log2FoldChange <= fc),])
  }
  else{
    stop("Please select 'up' or 'down' to query over-represented pathways among
              up- andn down-regulated genes, respectively.")
  }
  
  # KEGG module over-representation analysis
  # uses enricher
  mkk <- enrichMKEGG(gene = sig.genes,
                     organism = organism,
                     pvalueCutoff = p.cutoff,
                     pAdjustMethod = p.correction,
                     minGSSize = 1,
                     maxGSSize = 1e9)
  
  # KEGG module gene set enrichment analysis
  # rank genes by log2 fold change from highest to lowest
  genes.og <- genes.og[order(genes.og$log2FoldChange, decreasing = T),]
  allgenes.ranked <- genes.og$log2FoldChange
  names(allgenes.ranked) <- rownames(genes.og)
  print(allgenes.ranked)[1:5]
  
  # use GSEA function
  mkk2 <- gseMKEGG(geneList = allgenes.ranked,
                   organism = organism,
                   minGSSize = 1,
                   maxGSSize = 1e9,
                   eps = 0,
                   pvalueCutoff = p.cutoff,
                   pAdjustMethod = p.correction,
                   verbose = F)
  
  # get significant results with adjusted p<0.05
  if (!(is.null(mkk))){
  mkk.final <- mkk@result[which(mkk@result$p.adjust < 0.05),]}
  else{
    mkk.final <- data.frame()}
  if (!(is.null(mkk2))){
  mkk2.final <- mkk2@result[which(mkk2@result$p.adjust < 0.05),]}
  else{
    mkk2.final <- data.frame()}
  
  # add treatment column to final data frame
  mkk.final$treatment <- rep(treatment, nrow(mkk.final))
  mkk2.final$treatment <- rep(treatment, nrow(mkk2.final))
  
  # return final data frames
  return(list("module.overrep" = mkk.final,
              "module.gsea" = mkk2.final,
              "OGmodule.overrep.og" = mkk@result,
              "OGmodule.gsea.og" = mkk2@result))
}

###################### load/prepare kegg annotations

# clusterProfiler has LT2, so nothing else needed!
search_kegg_organism("stm", by="kegg_code")


###################### kegg over-representation analysis

# load annotated DESeq data
deseq.ultra <- read.delim(file = "11_deseq_final/annot_deseq_lt2.tsv",
                          header = T, sep = "\t",
                          stringsAsFactors = F,
                          check.names = F)
head(deseq.ultra)
rownames(deseq.ultra) <- unlist(lapply(strsplit(x = rownames(deseq.ultra), split = "_"), "[[", 1))
head(deseq.ultra)

# remove rows (genes) where log2 fold change is NA
deseq.ultra <- deseq.ultra[!(is.na(deseq.ultra$log2FoldChange)),]

# perform kegg pathway enrichment analysis (upregulated genes)
deseq.ultra.up.pathway <- run.kegg.pathway(genes.og = deseq.ultra,
                                           updown = "up", fc = 1.0,
                                           organism = "stm",
                                           p.cutoff = 0.05,
                                           p.correction = "fdr",
                                           treatment = "group")

dim(deseq.ultra.up.pathway$pathway.overrep)
dim(deseq.ultra.up.pathway$pathway.gsea)

# perform kegg pathway enrichment analysis (downregulated genes)
deseq.ultra.down.pathway <- run.kegg.pathway(genes.og = deseq.ultra,
                                             organism = "stm",
                                             updown = "down", fc = -1.0,
                                             p.cutoff = 0.05,
                                             p.correction = "fdr",
                                             treatment = "group")

dim(deseq.ultra.down.pathway$pathway.overrep)
dim(deseq.ultra.down.pathway$pathway.gsea)

# perform kegg module enrichment analysis (upregulated genes)
deseq.ultra.up.module <- run.kegg.module(genes.og = deseq.ultra,
                                         organism = "stm",
                                         updown = "up", fc = 1.0,
                                         p.cutoff = 0.05,
                                         p.correction = "fdr",
                                         treatment = "group")

dim(deseq.ultra.up.module$module.overrep)
dim(deseq.ultra.up.module$module.gsea)

# perform kegg module enrichment analysis (downregulated genes)
deseq.ultra.down.module <- run.kegg.module(genes.og = deseq.ultra,
                                           organism = "stm",
                                           updown = "down", fc = -1.0,
                                           p.cutoff = 0.05,
                                           p.correction = "fdr",
                                           treatment = "group")

dim(deseq.ultra.down.module$module.overrep)
dim(deseq.ultra.down.module$module.gsea)

###################### create and save final TSV files

# make data frame of enriched kegg pathways for upregulated genes
final.pathway.up <- rbind(deseq.ultra.up.pathway$pathway.overrep)

# save TSV file
write.table(x = final.pathway.up, file = "14_KEGG_enrichment/salmonella_kegg_pathway_up.tsv",
            append = F, quote = F, sep = "\t",
            row.names = T, col.names = T)

# make data frame of enriched kegg pathways for downregulated genes
final.pathway.down <- rbind(deseq.ultra.down.pathway$pathway.overrep)

# save TSV file
write.table(x = final.pathway.down, file = "14_KEGG_enrichment/salmonella_kegg_pathway_down.tsv",
            append = F, quote = F, sep = "\t",
            row.names = T, col.names = T)

# create final GSEA data frame for kegg pathways
final.pathway.gsea <- rbind(deseq.ultra.up.pathway$pathway.gsea)

# save TSV file
write.table(x = final.pathway.gsea, file = "14_KEGG_enrichment/salmonella_kegg_pathway_gsea.tsv",
            append = F, quote = F, sep = "\t",
            row.names = T, col.names = T)

# create final data frame with kegg module enrichment results for upregulated genes
final.module.up <- rbind(deseq.ultra.up.module$module.overrep)

# save TSV file
write.table(x = final.module.up, file = "14_KEGG_enrichment/salmonella_kegg_module_up.tsv",
            append = F, quote = F, sep = "\t",
            row.names = T, col.names = T)

# create final data frame with kegg module enrichment results for downregulated genes
final.module.down <- rbind(deseq.ultra.down.module$module.overrep)

# save TSV file
write.table(x = final.module.down, file = "14_KEGG_enrichment/salmonella_kegg_module_down.tsv",
            append = F, quote = F, sep = "\t",
            row.names = T, col.names = T)

# create final GSEA data frame for kegg modules
final.module.gsea <- rbind(deseq.ultra.up.module$module.gsea)

# save TSV file
write.table(x = final.module.gsea, file = "14_KEGG_enrichment/salmonella_kegg_module_gsea.tsv",
            append = F, quote = F, sep = "\t",
            row.names = T, col.names = T)

###################### create and save plots

# for kegg pathways, take -log10 of the raw p-values
final.pathway.up$log10p <- -log10(final.pathway.up$pvalue)
final.pathway.down$log10p <- log10(final.pathway.down$pvalue)

# merge overrepresented pathways
final.pathway.overrep <- rbind(final.pathway.up, final.pathway.down)
final.pathway.overrep$direction <- c(rep("Up-regulated", nrow(final.pathway.up)),
                                     rep("Down-regulated", nrow(final.pathway.down)))
head(final.pathway.overrep)

# create KEGG_term column
final.pathway.overrep$KEGG_term <- final.pathway.overrep$Description

# overview of results
table(final.pathway.overrep$direction)
summary(final.pathway.overrep$log10p)

# Reduce label text
head(final.pathway.overrep)
final.pathway.overrep$KEGG_term <- gsub(pattern = " - Salmonella enterica subsp. enterica serovar Typhimurium LT2",
                                        replacement = "", x = final.pathway.overrep$KEGG_term)
head(final.pathway.overrep)

# make bar chart of enriched kegg pathways
pdf(file = "14_KEGG_enrichment/plot_pathway_overrep.pdf",width = 8.5, height = 6)
plot.pathway.overrep <- ggplot(final.pathway.overrep[which(final.pathway.overrep$p.adjust<0.05),], aes(x=reorder(KEGG_term, log10p), y=log10p)) +#, fill = direction)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("KEGG Pathway") +
  ylab("(-)log10 P-Value") +
  #scale_y_continuous(limits = c(-5, 8), breaks = seq(-5, 8, by = 1)) +
  #scale_fill_manual(values = c("#228833", "#CCBB44", "#AA3377")) +
  theme_bw(base_size=12) +
  #guides(colour=guide_legend(override.aes=list(size=2.5))) +
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
plot.pathway.overrep
dev.off()

# for kegg modules, take -log10 of the raw p-values
final.module.up$log10p <- -log10(final.module.up$pvalue)
final.module.down$log10p <- log10(final.module.down$pvalue)

# merge overrepresented modules
final.module.overrep <- rbind(final.module.up, final.module.down)
final.module.overrep$direction <- c(rep("Up-regulated", nrow(final.module.up)),
                                    rep("Down-regulated", nrow(final.module.down)))
head(final.module.overrep)

# create KEGG_term column
final.module.overrep$KEGG_term <- final.module.overrep$Description

# overview of results
table(final.module.overrep$direction)
summary(final.module.overrep$log10p)

# make bar chart of enriched kegg modules
pdf(file = "14_KEGG_enrichment/plot_module_overrep.pdf",width = 8.5, height = 6)
plot.module.overrep <- ggplot(final.module.overrep[which(final.module.overrep$p.adjust<0.05),], aes(x=reorder(KEGG_term, log10p), y=log10p)) + #, fill = treatment)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("KEGG Module") +
  ylab("(-)log10 P-Value") +
  #scale_y_continuous(limits = c(-4, 5), breaks = seq(-4, 5, by = 1)) +
  #scale_fill_manual(values = c("#228833", "#CCBB44", "#AA3377")) +
  theme_bw(base_size=12) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
plot.module.overrep
dev.off()

# add KEGG_term column to final.pathway.gsea
final.pathway.gsea$KEGG_term <- final.pathway.gsea$Description

# Reduce label text
head(final.pathway.gsea)
final.pathway.gsea$KEGG_term <- gsub(pattern = " - Salmonella enterica subsp. enterica serovar Typhimurium LT2",
                                        replacement = "", x = final.pathway.gsea$KEGG_term)
head(final.pathway.gsea)


# plot GSEA results
pdf(file = "14_KEGG_enrichment/plot_pathway_gsea.pdf",width = 8.5, height = 6)
plot.pathway.gsea <- ggplot(final.pathway.gsea[which(final.pathway.gsea$p.adjust<0.05),], aes(x=reorder(KEGG_term, NES), y=NES)) + #, fill = treatment)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("KEGG Pathway") +
  ylab("NES") +
  #scale_y_continuous(limits = c(-2, 3), breaks = seq(-2, 3, by = 1)) +
  #scale_fill_manual(values = c("#228833", "#CCBB44", "#AA3377")) +
  theme_bw(base_size=12) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
plot.pathway.gsea
dev.off()

# add KEGG_term column to final.module.gsea
final.module.gsea$KEGG_term <- final.module.gsea$Description

# plot GSEA results
pdf(file = "14_KEGG_enrichment/plot_module_gsea.pdf",width = 8.5, height = 6)
plot.module.gsea <- ggplot(final.module.gsea[which(final.module.gsea$p.adjust<0.05),], aes(x=reorder(KEGG_term, NES), y=NES)) + #, fill = treatment)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("KEGG Module") +
  ylab("NES") +
  #scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1)) +
  #scale_fill_manual(values = c("#228833", "#CCBB44", "#AA3377")) +
  theme_bw(base_size=12) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  theme(panel.border = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()
plot.module.gsea
dev.off()
