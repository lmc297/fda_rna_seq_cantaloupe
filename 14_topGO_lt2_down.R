library(topGO)
library(ggplot2)
library(viridis)
library(gridExtra)

####################################### set seed

set.seed(42)

####################################### define GO enrichment function

# define run.fisher function, which will be used to perform GO enrichment 
run.fisher <- function(obj, algorithm, topNodes, pval, correct=NULL){
  
  # run Fisher's exact test using specified algorithm
  results.fisher <- runTest(obj, 
                            algorithm=algorithm,
                            statistic="fisher")
  
  # get the number of tests performed
  ntest <- length(results.fisher@score)
  
  # correct for multiple comparisons, if desired
  # you don't need to do this with topGO
  if (!(is.null(correct))){
    results.fisher@score <- p.adjust(p = results.fisher@score, method = correct)
  }
  
  # perform enrichment using specified topNodes
  # topNodes == number of top GO terms to be included in the table
  goEnrichment.full <- GenTable(obj, 
                                Test=results.fisher, 
                                orderBy="Test", 
                                topNodes=topNodes)
  
  # note: you might need to convert values with very low p-values, which topGo does not treat as numeric
  print(goEnrichment.full[grepl(pattern = "<", x = goEnrichment.full$Test),])
  goEnrichment.full$Test[grepl(pattern = "<", x = goEnrichment.full$Test)] <- "1e-30"
  goEnrichment.full$Test <- as.numeric(goEnrichment.full$Test)
  
  # get significant p-values
  goEnrichment <- goEnrichment.full
  goEnrichment$Test <- as.numeric(goEnrichment$Test)
  goEnrichment <- goEnrichment[goEnrichment$Test<pval,]
  
  # make a pretty picture
  goEnrichment <- goEnrichment[,c("GO.ID","Term","Test")]
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
  
  p <- ggplot(goEnrichment, aes(x=Term, y=-log10(Test))) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab("Gene Ontology (GO) Term") +
    ylab("-log10(P-Value)") +
    scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$Test)), by = 2), 1)) +
    theme_bw(base_size=24) +
    theme(
      legend.position='none',
      legend.background=element_rect(),
      axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.10),
      axis.text.y=element_text(angle=0, size=12, face="bold", vjust=0.5),
      axis.title=element_text(size=24, face="bold"),
      legend.key=element_blank(),     
      legend.key.size=unit(1, "cm"),      
      legend.text=element_text(size=18),  
      title=element_text(size=18)) +
    guides(colour=guide_legend(override.aes=list(size=2.5))) +
    coord_flip()
  
  # return multiple elements
  return(list("df" = goEnrichment.full, "sig.df" = goEnrichment, "plot" = p, "ntest" = ntest))
  
  
}

####################################### load/clean input data from DESeq

# load annotated DESeq results
crf.og <- read.delim(file = "11_deseq_final/annot_deseq_lt2.tsv", 
                     header = T, sep = "\t", 
                     stringsAsFactors = F, 
                     check.names = F, row.names = 1)
head(crf.og)
dim(crf.og)

# subset annotated DESeq results
# take locus name, log2 fold change, raw p-value, adjusted p-value
crf <- data.frame(unlist(lapply(strsplit(x = rownames(crf.og), split = "_"), "[[", 1)),
                  crf.og$log2FoldChange, crf.og$pvalue, crf.og$padj)
head(crf)
colnames(crf) <- c("locus", "fc", "pvalue", "padj")
head(crf)

# remove duplicated rows
crf <- crf[!duplicated(crf), ]
dim(crf)

# read gene2GO mappings, just for fun
geneID2GO <- readMappings(file = "12_GO_terms/eggnog.tsv")
str(head(geneID2GO))
geneID2GO[["locus"]] <- NULL
str(head(geneID2GO))
length(geneID2GO)

# remove unused
geneID2GO <- geneID2GO[unlist(lapply(X = geneID2GO, FUN = function(x) length(x)>0))]
crf <- crf[which(crf$locus%in%names(geneID2GO)),]
names(geneID2GO)[!(names(geneID2GO)%in%crf$locus)]

# sanity check (see if all genes are accounted for)
table(names(geneID2GO)%in%crf$locus)
table(crf$locus%in%names(geneID2GO))

####################################### get gene rankings

# rank genes by raw p-value 
crf.ranked <- crf[order(crf$pvalue, decreasing = F),]
head(crf.ranked)
crf.ranked[1:30,]

# save gene names 
crf.ranked.names <- crf.ranked$locus

# set up ranked genes for fisher's exact test
# we consider genes to be significant if adjusted p-value < 0.05 AND  fold change <=-1
crf.ranked <- as.factor(as.integer(ifelse(test = (crf.ranked$padj<0.05 & crf.ranked$fc<=-1.0),
                                          yes = 1, no = 0)))
names(crf.ranked) <- crf.ranked.names
head(crf.ranked)
table(crf.ranked)

####################################### test biological process

# create new topGO object with biological process
bgc.bp <- new("topGOdata", ontology="BP", allGenes=crf.ranked,
              annotationFun=annFUN.file,
              file = "12_GO_terms/eggnog.tsv",
              nodeSize=3)

# run topGO enrichment with weight01 algorithm and all nodes
bgc.bp.result <- run.fisher(obj = bgc.bp, algorithm = "weight01",
                            topNodes = 1699, pval = 0.05)


dim(bgc.bp.result$df)
dim(bgc.bp.result$sig.df)
bgc.bp.result$ntest

# plot results
pdf(file = "12_GO_terms/go_lt2_bp_down.pdf", width = 8.5, height = 7)
bgc.bp.result$plot
dev.off()

####################################### test molecular function

# create new topGO object with molecular function
bgc.mf <- new("topGOdata", ontology="MF", allGenes=crf.ranked,
              annotationFun=annFUN.file,
              file = "12_GO_terms/eggnog.tsv",
              nodeSize=3)

# run topGO enrichment with weight01 algorithm and all nodes
bgc.mf.result <- run.fisher(obj = bgc.mf, algorithm = "weight01",
                            topNodes = 794, pval = 0.05)

dim(bgc.mf.result$df)
dim(bgc.mf.result$sig.df)
bgc.mf.result$ntest

# plot results
pdf(file = "12_GO_terms/go_lt2_mf_down.pdf", width = 8.5, height = 7)
bgc.mf.result$plot
dev.off()

####################################### test cellular component

# create new topGO object with cellular component
bgc.cc <- new("topGOdata", ontology="CC", allGenes=crf.ranked,
              annotationFun=annFUN.file,
              file = "12_GO_terms/eggnog.tsv",
              nodeSize=3)

# run topGO enrichment with weight01 algorithm and all nodes
bgc.cc.result <- run.fisher(obj = bgc.cc, algorithm = "weight01",
                            topNodes = 159, pval = 0.05)

dim(bgc.cc.result$df)
dim(bgc.cc.result$sig.df)
bgc.cc.result$ntest

# plot results
pdf(file = "12_GO_terms/go_lt2_cc_down.pdf", width = 8.5, height = 6)
bgc.cc.result$plot
dev.off()

####################################### create final data frame

# create final data frame with all results
final.df <- rbind(bgc.bp.result$df,
                  bgc.mf.result$df,
                  bgc.cc.result$df)
head(final.df)

# add ontology column
final.df$Ontology <- c(rep("BP", nrow(bgc.bp.result$df)),
                       rep("MF", nrow(bgc.mf.result$df)),
                       rep("CC", nrow(bgc.cc.result$df)))
table(final.df$Ontology)

# save results as TSV
write.table(x = final.df, file = "12_GO_terms/go_lt2_down.tsv",
            append = F, quote = F,
            sep = "\t", row.names = F, col.names = T)
