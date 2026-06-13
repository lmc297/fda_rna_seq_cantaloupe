library(ggplot2)
library(stringr)
library(reshape2)
library(ape)
library(khroma)
library(Rsubread)

# inspired by BactSeq: https://github.com/adamd3/BactSeq/blob/main/bin/count_reads.R

set.seed(12345)

############# Read in reference genome annotation GFF

# read reference GFF file from NCBI
ref_annot <- ape::read.gff(file = "reference_genome_lt2/salmonella_GCF_000006945.2.gff",
                           na.strings = c(".", "?"), GFF3 = T)

# see what we're working with
head(ref_annot)
table(ref_annot$type)

# keep features with type == "gene"
ref_annot <- subset(ref_annot, type == "gene")

# get gene attributes 
gene_attr <- stringr::str_split(ref_annot$attributes, ";")

locus_tags <- unlist(lapply(gene_attr, function(x) {
  x[grepl("locus_tag", x)]
}))

gene_biotypes <- unlist(lapply(gene_attr, function(x) {
  x[grepl("gene_biotype", x)]
}))

common_gene_names <- unlist(lapply(gene_attr, function(x) {
  x <- x[grepl("gene=", x)]
  x[identical(x, character(0))] <- ""
  x
}))

# get gene lengths
gene_lengths <- (ref_annot$end - ref_annot$start) + 1

# make ref_gene_df, a data frame of gene attributes + lengths
ref_gene_df <- data.frame(
  locus_tag = locus_tags,
  biotype = gene_biotypes,
  gene_name = common_gene_names,
  gene_length = gene_lengths
)
ref_gene_df$locus_tag <- gsub("locus_tag=", "", ref_gene_df$locus_tag)
ref_gene_df$biotype <- gsub("gene_biotype=", "", ref_gene_df$biotype)
ref_gene_df$gene_name <- gsub("gene=", "", ref_gene_df$gene_name)

# take a look at ref_gene_df
head(ref_gene_df)
table(ref_gene_df$biotype)

# save ref_gene_df as a TSV file
write.table(x = ref_gene_df, file = "reference_genome_lt2/ref_gene_df.tsv",
            col.names = T, row.names = F,
            sep = "\t", quote = F)

############# Make read count tables

# collect BAM files
bamfiles <- Sys.glob(paths = "bam_pass2/*.bam")
bamfiles

# count reads using featureCounts
gene_counts <- Rsubread::featureCounts(
  files = bamfiles,
  annot.ext = "reference_genome_lt2/salmonella_GCF_000006945.2.gff",
  isGTFAnnotationFile = T,
  GTF.featureType = "gene",
  GTF.attrType = "locus_tag",
  nthreads = 1,
  countMultiMappingReads = F,
  isPairedEnd = T,
  strandSpecific = 2
)

# create stats_mat, a data frame of stats produced by featureCounts
stats_mat <- as.data.frame(gene_counts$stat)
colnames(stats_mat) <- gsub(pattern = "_mapped.bam", replacement = "", x = colnames(stats_mat))
head(stats_mat)

# create counts_mat, a data frame of read counts produced by featureCounts
counts_mat <- as.data.frame(gene_counts$counts)
head(counts_mat)
colnames(counts_mat) <- gsub(pattern = "_mapped.bam", replacement = "", x = colnames(counts_mat))
head(counts_mat)

# sanity check: make sure ref_gene_df and counts_mat can be merged
table(rownames(counts_mat)%in%ref_gene_df$locus_tag)
table(ref_gene_df$locus_tag%in%rownames(counts_mat))
table(rownames(counts_mat)==ref_gene_df$locus_tag)
table(ref_gene_df$locus_tag==rownames(counts_mat))

# merge ref_gene_df and counts_mat into counts_mat
counts_mat <- cbind(ref_gene_df, counts_mat)
head(counts_mat)

# save counts_mat as TSV file
write.table(x = counts_mat, file = "9_count_reads/gene_counts.tsv",
            col.names = T, row.names = F,
            sep = "\t", quote = F)

# subset counts_mat to get protein-coding genes only
gene_counts_pc <- counts_mat[counts_mat$biotype == "protein_coding", ]
head(gene_counts_pc)

# save counts for protein-coding genes as TSV file
write.table(x = gene_counts_pc, file = "9_count_reads/gene_counts_pc.tsv",
            col.names = T, row.names = F,
            sep = "\t", quote = F)

############# Summarize read counts per biotype

# get a vector of gene biotypes for the reference genome
all_biotypes <- unique(counts_mat$biotype)
all_biotypes

# create biotype_counts, a data frame of read counts per biotype
biotype_counts <- data.frame(do.call(
  cbind,
  lapply(all_biotypes, function(biotype) {
    colSums(gene_counts$counts[counts_mat$biotype == biotype, , drop = FALSE])
  })
))
colnames(biotype_counts) <- all_biotypes
rownames(biotype_counts) <- gsub(pattern = "_mapped.bam", replacement = "", x = rownames(biotype_counts))
head(biotype_counts)

# make counts_summary, a data frame with sample metadata
counts_summary <- data.frame(rownames(biotype_counts),
                             unlist(lapply(strsplit(x = rownames(biotype_counts), split = "_"), "[[", 4)),
                             unlist(lapply(strsplit(x = rownames(biotype_counts), split = "_"), "[[", 2)))
colnames(counts_summary) <- c("sample", "group", "rep_no")
counts_summary

# sanity check: make sure counts_summary and biotype_counts can be merged
table(counts_summary$sample==rownames(biotype_counts))
table(rownames(biotype_counts)==counts_summary$sample)

# merge counts_summary and biotype_counts into counts_summary
counts_summary <- cbind(counts_summary, biotype_counts)
counts_summary

# collect files with total # of mapped reads per sample
countfiles <- Sys.glob(paths = "bam_pass1/*.counts")
countfiles

# create merged_total_counts, a data frame with number of mapped reads per sample
fvec <- c()
cvec <- c()
for (i in 1:length(countfiles)){
  fname <- countfiles[i]
  fname <- gsub(pattern = "bam_pass1/", replacement = "", x = fname)
  fname <- gsub(pattern = ".counts", replacement = "", x = fname)
  fvec <- c(fvec, fname)
  cf <- read.delim(file = countfiles[i],
                   header = F, sep = "\t",
                   stringsAsFactors = F,
                   check.names = F)
  cvec <- c(cvec, cf$V1)
}
merged_total_counts <- data.frame(fvec, cvec)
colnames(merged_total_counts) <- c("sample", "mapped")
merged_total_counts

# sanity check: make sure counts_summary and merged_total_counts can be merged
table(counts_summary$sample%in%merged_total_counts$sample)
table(merged_total_counts$sample%in%counts_summary$sample)
table(counts_summary$sample==merged_total_counts$sample)
table(merged_total_counts$sample==counts_summary$sample)

# reorder merged_total_counts to match counts_summary and sanity check
merged_total_counts <- merged_total_counts[match(counts_summary$sample, merged_total_counts$sample),]
table(counts_summary$sample==merged_total_counts$sample)
table(merged_total_counts$sample==counts_summary$sample)

# merge counts_summary and merged_total_counts into counts_summary
counts_summary$mapped <- merged_total_counts$mapped
counts_summary

# create missing column (total number of mapped reads, minus reads assigned to any biotype)
counts_summary$missing <- counts_summary$mapped - rowSums(counts_summary[,all_biotypes])
counts_summary
counts_summary$missing[which(counts_summary$missing<0)] <- 0
counts_summary

# create nonprot column (total number of mapped reads, minus reads assigned to protein-coding genes)
counts_summary$nonprot <- counts_summary$mapped - counts_summary$protein_coding
counts_summary

# save counts_summary as a TSV file
write.table(x = counts_summary, file = "9_count_reads/counts_summary.tsv",
            col.names = T, row.names = F,
            sep = "\t", quote = F)

############# Plot read counts for protein-coding genes vs everything else (nonprot column)

# create counts_melt, a data frame with read counts for protein-coding genes vs everything else (nonprot column)
counts_melt <- reshape2::melt(
  counts_summary,
  id.vars = c("sample"),
  measure.vars = c(
    "protein_coding",
    "nonprot"
  )
)

counts_melt

# plot counts_melt (total read counts assigned to protein-coding genes vs everything else)
p1 <- ggplot(
  counts_melt,
  aes(x = sample, colour = variable, fill = variable, y = value)
) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  coord_flip() +
  xlab("Sample") +
  ylab("# read pairs") + 
  theme_bw() + 
  scale_fill_manual(values = c(protein_coding = "turquoise3", nonprot = "deeppink3")) +
  scale_color_manual(values = c(protein_coding = "turquoise3", nonprot = "deeppink3"))
p1

# create propCols, a data frame of the proportion of reads assigned to protein-coding genes vs everything else (nonprot column)
propCols <- counts_summary[c(
  "protein_coding",
  "nonprot"
)] / counts_summary$mapped # divide by total number of mapped reads for each sample

# sanity check: each row should sum to 1
rowSums(propCols) 

# add sample name column to propCols
propCols$sample <- counts_summary$sample

# create prop_melt, a data frame with proportion of reads assigned to protein-coding genes vs everything else (nonprot column)
prop_melt <- reshape2::melt(
  propCols,
  id.vars = c("sample"),
  measure.vars = c(
    "protein_coding",
    "nonprot"
  )
)

# plot prop_melt (proportion of reads assigned to protein-coding genes vs everything else)
p2 <- ggplot(
  prop_melt,
  aes(x = sample, colour = variable, fill = variable, y = value)
) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  xlab("Sample") +
  ylab("Proportion of reads") + 
  theme_bw() + 
  scale_fill_manual(values = c(protein_coding = "turquoise3", nonprot = "deeppink3")) +
  scale_color_manual(values = c(protein_coding = "turquoise3", nonprot = "deeppink3"))
p2

# save plots p1 and p2 as PDF files
pdf(file = "9_count_reads/bar_counts_total.pdf", width = 8.5, height = 6)
p1
dev.off()

pdf(file = "9_count_reads/bar_counts_proportion.pdf", width = 8.5, height = 6)
p2
dev.off()

############# Plot read counts for all biotypes vs everything else (missing column)

# create counts_melt2, a data frame with read counts for all biotypes vs everything else (missing column)
counts_melt2 <- reshape2::melt(
  counts_summary,
  id.vars = c("sample"),
  measure.vars = c(
    all_biotypes,
    "missing"
  )
)

counts_melt2
table(counts_melt2$variable)

# plot counts_melt2 (total read counts assigned to all biotypes vs everything else)
muted <- color("muted")
p3 <- ggplot(
  counts_melt2,
  aes(x = sample, colour = variable, fill = variable, y = value)
) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  coord_flip() +
  xlab("Sample") +
  ylab("# read pairs") + 
  theme_bw() + 
  scale_fill_manual(values = muted(n = 9)) +
  scale_color_manual(values = muted(n = 9))
p3

# create propCols.expanded, a data frame of the proportion of reads assigned to each biotype vs everything else (missing column)
propCols.expanded <- counts_summary[c(
  all_biotypes, 
  "missing"
)] / counts_summary$mapped # divide by total number of mapped reads for each sample

# sanity check: each row should sum to 1
rowSums(propCols.expanded)

# add sample name column to propCols.expanded
propCols.expanded$sample <- counts_summary$sample

# create prop_melt.expanded, a data frame with proportion of reads assigned to all biotypes vs everything else (missing column)
prop_melt.expanded <- reshape2::melt(
  propCols.expanded,
  id.vars = c("sample"),
  measure.vars = c(
    all_biotypes,
    "missing"
  )
)

# plot prop_melt.expanded (proportion of reads assigned to all biotypes vs everything else)
p4 <- ggplot(
  prop_melt.expanded,
  aes(x = sample, colour = variable, fill = variable, y = value)
) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  xlab("Sample") +
  ylab("Proportion of reads") + 
  theme_bw() + 
  scale_fill_manual(values = muted(n = 9)) +
  scale_color_manual(values = muted(n = 9))
p4

# save plots p3 and p4 as PDF files
pdf(file = "9_count_reads/bar_counts_total_expanded.pdf", width = 8.5, height = 6)
p3
dev.off()

pdf(file = "9_count_reads/bar_counts_proportion_expanded.pdf", width = 8.5, height = 6)
p4
dev.off()

############# Plot featureCounts statistics

# transpose stats_mat so that samples are in rows and counts are in columns
stats_mat = setNames(data.frame(t(stats_mat[,-1])), stats_mat[,1])
stats_mat

# remove columns that only contain zeros
stats_mat <- stats_mat[,which(colSums(stats_mat)>0)]
stats_mat

# add a column containing sample names
stats_mat$sample <- rownames(stats_mat)
stats_mat

# create stats_melt, a data frame with read counts per featureCounts stat
stats_melt <- reshape2::melt(
  stats_mat,
  id.vars = c("sample"),
  measure.vars = c(
    colnames(stats_mat[1:(ncol(stats_mat)-1)])
  )
)

stats_melt
table(stats_melt$variable)

# plot stats_melt (total read counts assigned to each featureCounts stat)
vibrant <- color("vibrant")
p5 <- ggplot(
  stats_melt,
  aes(x = sample, colour = variable, fill = variable, y = value)
) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  coord_flip() +
  xlab("Sample") +
  ylab("# read pairs") + 
  theme_bw() + 
  scale_fill_manual(values = vibrant(n = 3)) +
  scale_color_manual(values = vibrant(n = 3))
p5

# create statCols, a data frame of the proportion of reads assigned to each featureCounts stat
# divide by total number of featureCounts stat assignments per sample
# note that this is NOT the same thing as total number of mapped reads per sample
statCols <- stats_mat[c(
  colnames(stats_mat)[1:(ncol(stats_mat)-1)]
)] / rowSums(stats_mat[,1:(ncol(stats_mat)-1)]) # divide by total number of featureCounts stat assignments

# sanity check: each row should sum to 1
rowSums(statCols)

# add column name with sample
statCols$sample <- stats_mat$sample
statCols

# create propStat_melt, a data frame with proportion of reads assigned to each featureCounts stat
propStat_melt <- reshape2::melt(
  statCols,
  id.vars = c("sample"),
  measure.vars = c(
    colnames(statCols[1:(ncol(statCols)-1)])
  )
)
propStat_melt

# plot propStat_melt (proportion of reads assigned to each featureCounts stat)
p6 <- ggplot(
  propStat_melt,
  aes(x = sample, colour = variable, fill = variable, y = value)
) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  xlab("Sample") +
  ylab("Proportion of reads") + 
  theme_bw() + 
  scale_fill_manual(values = vibrant(n = 3)) +
  scale_color_manual(values = vibrant(n = 3))
p6

# save plots p5 and p6 as PDF files
pdf(file = "9_count_reads/bar_stats_total.pdf", width = 8.5, height = 6)
p5
dev.off()

pdf(file = "9_count_reads/bar_stats_proportion.pdf", width = 8.5, height = 6)
p6
dev.off()
