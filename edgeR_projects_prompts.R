#!/usr/bin/env Rscript

# R script to perform statistical analysis of gene count tables using edgeR exact tests
# usage: Rscript edgeR_prompts.r mergedCountsFile
# usage Ex: Rscript edgeR_prompts.r /YOUR/PATH/Adipocyte/counted/merged_counts_formatted.txt
# usage Ex: Rscript edgeR_prompts.r /YOUR/PATH/Jurkat/counted/merged_counts_formatted.txt


#
## Setup Stage
#

# set working directory
setwd("C:/Users/bhava/OneDrive/Desktop/Scripts")

# install packages, this should only need to be done once
#Install edgeR, this should only need to be done once
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# #Install packages, this should only need to be done once
# BiocManager::install("edgeR", force = TRUE)
# BiocManager::install("biomaRt", force = TRUE)
# install.packages("msigdbr")
# install.packages("ggplot2")
# install.packages("tibble")
# install.packages("dplyr")
# install.packages("pheatmap")
# install.packages("ggplotify")

# load libraries
library(edgeR)
library(biomaRt)
library(msigdbr)
library(ggplot2)
library(tibble)
library(dplyr)
library(pheatmap)
library(ggplotify)


# connect to a an ensembl website mart by specifying a BioMart and dataset parameters for Mus_musculus.GRCm39
#Mus_musculus.GRCm39
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
#Homo_sapiens.GRCh38
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")


# retrieve gene attributes
genes_att <- getBM(attributes=c('ensembl_gene_id',
                                'external_gene_name',
                                #'uniprot_gn_id',
                                'description'), 
                   mart = ensembl)


# retrieve all hallmark gene sets
#h_gene_sets = msigdbr(species = "mouse", category = "H")
h_gene_sets = msigdbr(species = "human", category = "H")


# retrieve the apoptosis gene set
h_apoptosis <- h_gene_sets %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS")


# retrieve the TNF alpha signaling gene set
h_tnfa <- h_gene_sets %>%
  dplyr::filter(gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB")


# retrieve the inflammatory response gene set
h_inflammatory <- h_gene_sets %>%
  dplyr::filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE")


# retrieve the interferon alpha response gene set
h_interferon <- h_gene_sets %>%
  dplyr::filter(gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE")


# import gene count data
#inputTable <- read.table(file="Adipocyte_merged_counts_formatted.txt", header = TRUE, sep = "\t", row.names="gene")
inputTable <- read.table(file="Jurkat_merged_counts_formatted.txt", header = TRUE, sep = "\t", row.names="gene")


# subset the input counts
#DE Adipocyte
#subsetTable <- inputTable[ , -which(names(inputTable) %in% c("S00", "S01", "S06", "S10"))]
#DE Jurkat
subsetTable <- inputTable[ , -which(names(inputTable) %in% c("S00", "S04", "S07", "S09"))]


# trim the data table
countsTable <- head(subsetTable, - 5)


# set number of samples
numSamples <- 9


# add grouping factor
group <- factor(c(rep("370mW/cm2",3), rep("444mW/cm2",3), rep("CTL",3)))


# create DGE list object
list <- DGEList(counts=countsTable,group=group)


#
## Analysis Prep Stage
#

# retrieve library sizes
libraries <- data.frame(
  samples = names(countsTable),
  sizes = list$samples$lib.size*1e-6)


# plot the library sizes before normalization
jpeg("plotBars_librarySizes.jpg")
ggplot(data = libraries, aes(x = samples, y = sizes)) + 
  geom_bar(stat="identity") +
  labs(x = "Sample", y="Library size (millions)")
dev.off()


# there is no purpose in analyzing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
                     

# keep only the genes that meet our criteria
list <- list[keep, , keep.lib.sizes=FALSE]


# calculate normalized factors
list <- calcNormFactors(list)
                        

# calculate cpm of normalized counts
normList <- cpm(list, normalized.lib.sizes=TRUE)


# write normalized cpm to file
write.table(normList, file="normalized_counts.csv", sep=",", row.names=TRUE)


# draw a MDS plot to show the relative similarities of the samples
# and to view batch and treatment effects after normalization
jpeg("plotMDS_afterNormalize.jpg")
plotMDS(list, col=rep(1:3, each=3))
dev.off()


# draw a heatmap of individual RNA-seq samples using moderated
# log-counts-per-million after normalization
logcpm <- cpm(list, log=TRUE)
jpeg("heatmap_afterNormalize.jpg")
heatmap(logcpm)
dev.off()


# produce a matrix of pseudo-counts using estimates of the common and tagwise dispersions
list <- estimateDisp(list)


# view dispersion estimates and biological coefficient of variation
jpeg("plotBCV.jpg")
plotBCV(list)
dev.off()


#
## DE Analysis Stage
#

# setup color blind safe color vectors for plotting
color_subset <- c("#0000FF", "#000000", "#FF0000")


# prepare a gene counts table
sample_counts <- data.frame(logcpm) %>% rownames_to_column(var="gene")


# create data frame with the experimental design layout
exp_factor <- data.frame(group)


# attach rownames to the data frame
rownames(exp_factor) <- colnames(logcpm)


##
### 370mW/cm2 vs CTL
##

# perform an exact test for 100mV vs CTL
tested <- exactTest(list, pair=c("CTL", "370mW/cm2"))


# create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table %>% rownames_to_column(var="gene")
                                                                                 

# write results table of DE genes to a CSV file
write.csv(resultsTbl, file="370mW_CTL_topTags.csv")


# plot log-fold change against log-counts per million, with DE genes highlighted
# the blue lines indicate 2-fold changes
jpeg("370mW_CTL_plotMD.jpg")
plotMD(tested)
abline(h=c(-1, 1), col="grey")
dev.off()


# identify significantly DE genes
resultsTbl$sigDE <- "NA"
resultsTbl$sigDE[resultsTbl$logFC > 1 & resultsTbl$FDR < 0.05] <- "UP"
resultsTbl$sigDE[resultsTbl$logFC < -1 & resultsTbl$FDR < 0.05] <- "DOWN"


# create volcano plot
jpeg("370mW_CTL_plotVolcano.jpg")
ggplot(data=resultsTbl, aes(x=logFC, y=-log10(FDR), color = sigDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = color_subset)
dev.off()


# remove the sigDE column
resultsTbl_clean <- resultsTbl[ , -which(names(resultsTbl) %in% c("sigDE"))]


# return all rows from resultsTbl_clean and genes_att
resultsTbl_att <- resultsTbl_clean %>% left_join(genes_att, by = c("gene" = "ensembl_gene_id"))


# sort the results by FDR
resultsTbl_ann <- resultsTbl_att[order(resultsTbl_att$FDR),]


# write the sorted DE genes to a CSV file
write.csv(resultsTbl_ann, file="370mW_CTL_annotated.csv")


# filter for significantly DE genes
resultsTbl_filt <- resultsTbl_ann[resultsTbl_ann$FDR < 0.05,]


# return all rows from resultsTbl_filt and sample_counts
resultsTbl_counts <- resultsTbl_filt %>% inner_join(sample_counts, by = "gene")


# return all rows from resultsTbl with a match in h_apoptosis and format for plotting
resultsTbl_apoptosis <- resultsTbl_counts %>% semi_join(h_apoptosis, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_apoptosis_filt <- resultsTbl_apoptosis[ , -which(names(resultsTbl_apoptosis) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_apoptosis_filt) <- resultsTbl_apoptosis$external_gene_name
resultsTbl_apoptosis_filt <- head(resultsTbl_apoptosis_filt, n = 20)


# return all rows from resultsTbl with a match in h_inflammatory and format for plotting
resultsTbl_inflammatory <- resultsTbl_counts %>% semi_join(h_inflammatory, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_inflammatory_filt <- resultsTbl_inflammatory[ , -which(names(resultsTbl_inflammatory) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_inflammatory_filt) <- resultsTbl_inflammatory$external_gene_name
resultsTbl_inflammatory_filt <- head(resultsTbl_inflammatory_filt, n = 20)


# return all rows from resultsTbl with a match in h_tnfa and format for plotting
resultsTbl_tnfa <- resultsTbl_counts %>% semi_join(h_tnfa, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_tnfa_filt <- resultsTbl_tnfa[ , -which(names(resultsTbl_tnfa) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_tnfa_filt) <- resultsTbl_tnfa$external_gene_name
resultsTbl_tnfa_filt <- head(resultsTbl_tnfa_filt, n = 20)


# return all rows from resultsTbl with a match in h_interferon and format for plotting
resultsTbl_interferon <- resultsTbl_counts %>% semi_join(h_interferon, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_interferon_filt <- resultsTbl_interferon[ , -which(names(resultsTbl_interferon) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_interferon_filt) <- resultsTbl_interferon$external_gene_name
resultsTbl_interferon_filt <- head(resultsTbl_interferon_filt, n = 20)


# create heatmap for resultsTbl_apoptosis_filt
as.ggplot(pheatmap(resultsTbl_apoptosis_filt, scale="row", annotation_col = exp_factor,
                   main="370mW/cm2 vs CTL Apoptosis Response"))
ggsave("370mW_CTL_heatmap_apoptosis.png", bg = "white")


# create heatmap for resultsTbl_inflammatory_filt
as.ggplot(pheatmap(resultsTbl_inflammatory_filt, scale="row", annotation_col = exp_factor,
                   main="370mW/cm2 vs CTL Inflammatory Response"))
ggsave("370mW_CTL_heatmap_inflammatory.png", bg = "white")


# create heatmap for resultsTbl_tnfa_filt
as.ggplot(pheatmap(resultsTbl_tnfa_filt, scale="row", annotation_col = exp_factor,
                   main="370mW/cm2 vs CTL TNF Alpha Response"))
ggsave("370mW_CTL_heatmap_tnfa.png", bg = "white")


# create heatmap for resultsTbl_interferon_filt
as.ggplot(pheatmap(resultsTbl_interferon_filt, scale="row", annotation_col = exp_factor,
                   main="370mW/cm2 vs CTL Interferon Alpha Response"))
ggsave("370mW_CTL_heatmap_interferon.png", bg = "white")



##
### 444mW/cm2 vs CTL
##

# perform an exact test for 180mV vs CTL
tested <- exactTest(list, pair=c("CTL", "444mW/cm2"))

# create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table %>% rownames_to_column(var="gene")


# write results table of DE genes to a CSV file
write.csv(resultsTbl, file="444mW_CTL_topTags.csv")


# plot log-fold change against log-counts per million, with DE genes highlighted
# the blue lines indicate 2-fold changes
jpeg("444mW_CTL_plotMD.jpg")
plotMD(tested)
abline(h=c(-1, 1), col="grey")
dev.off()


# identify significantly DE genes
resultsTbl$sigDE <- "NA"
resultsTbl$sigDE[resultsTbl$logFC > 1 & resultsTbl$FDR < 0.05] <- "UP"
resultsTbl$sigDE[resultsTbl$logFC < -1 & resultsTbl$FDR < 0.05] <- "DOWN"


# create volcano plot
jpeg("444mW_CTL_plotVolcano.jpg")
ggplot(data=resultsTbl, aes(x=logFC, y=-log10(FDR), color = sigDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = color_subset)
dev.off()


# remove the sigDE column
resultsTbl_clean <- resultsTbl[ , -which(names(resultsTbl) %in% c("sigDE"))]


# return all rows from resultsTbl_clean and genes_att
resultsTbl_att <- resultsTbl_clean %>% left_join(genes_att, by = c("gene" = "ensembl_gene_id"))


# sort the results by FDR
resultsTbl_ann <- resultsTbl_att[order(resultsTbl_att$FDR),]


# write the sorted DE genes to a CSV file
write.csv(resultsTbl_ann, file="444mW_CTL_annotated.csv")


# filter for significantly DE genes
resultsTbl_filt <- resultsTbl_ann[resultsTbl_ann$FDR < 0.05,]


# return all rows from resultsTbl_filt and sample_counts
resultsTbl_counts <- resultsTbl_filt %>% inner_join(sample_counts, by = "gene")


# return all rows from resultsTbl with a match in h_apoptosis and format for plotting
resultsTbl_apoptosis <- resultsTbl_counts %>% semi_join(h_apoptosis, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_apoptosis_filt <- resultsTbl_apoptosis[ , -which(names(resultsTbl_apoptosis) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_apoptosis_filt) <- resultsTbl_apoptosis$external_gene_name
resultsTbl_apoptosis_filt <- head(resultsTbl_apoptosis_filt, n = 20)


# return all rows from resultsTbl with a match in h_inflammatory and format for plotting
resultsTbl_inflammatory <- resultsTbl_counts %>% semi_join(h_inflammatory, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_inflammatory_filt <- resultsTbl_inflammatory[ , -which(names(resultsTbl_inflammatory) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_inflammatory_filt) <- resultsTbl_inflammatory$external_gene_name
resultsTbl_inflammatory_filt <- head(resultsTbl_inflammatory_filt, n = 20)


# return all rows from resultsTbl with a match in h_tnfa and format for plotting
resultsTbl_tnfa <- resultsTbl_counts %>% semi_join(h_tnfa, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_tnfa_filt <- resultsTbl_tnfa[ , -which(names(resultsTbl_tnfa) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_tnfa_filt) <- resultsTbl_tnfa$external_gene_name
resultsTbl_tnfa_filt <- head(resultsTbl_tnfa_filt, n = 20)


# return all rows from resultsTbl with a match in h_interferon and format for plotting
resultsTbl_interferon <- resultsTbl_counts %>% semi_join(h_interferon, by = c("external_gene_name" = "gene_symbol"), na_matches = c("na", "never"))
resultsTbl_interferon_filt <- resultsTbl_interferon[ , -which(names(resultsTbl_interferon) %in% c("gene", "logFC", "logCPM", "PValue", "FDR", "ensembl_gene_id", "external_gene_name", "description"))]
rownames(resultsTbl_interferon_filt) <- resultsTbl_interferon$external_gene_name
resultsTbl_interferon_filt <- head(resultsTbl_interferon_filt, n = 20)


# create heatmap for resultsTbl_apoptosis_filt
as.ggplot(pheatmap(resultsTbl_apoptosis_filt, scale="row", annotation_col = exp_factor,
                   main="444mW/cm2 vs CTL Apoptosis Response"))
ggsave("444mW_CTL_heatmap_apoptosis.png", bg = "white")


# create heatmap for resultsTbl_inflammatory_filt
as.ggplot(pheatmap(resultsTbl_inflammatory_filt, scale="row", annotation_col = exp_factor,
                   main="444mW/cm2 vs CTL Inflammatory Response"))
ggsave("444mW_CTL_heatmap_inflammatory.png", bg = "white")


# create heatmap for resultsTbl_tnfa_filt
as.ggplot(pheatmap(resultsTbl_tnfa_filt, scale="row", annotation_col = exp_factor,
                   main="444mW/cm2 vs CTL TNF Alpha Response"))
ggsave("444mW_CTL_heatmap_tnfa.png", bg = "white")


# create heatmap for resultsTbl_interferon_filt
as.ggplot(pheatmap(resultsTbl_interferon_filt, scale="row", annotation_col = exp_factor,
                   main="444mW/cm2 vs CTL Interferon Alpha Response"))
ggsave("444mW_CTL_heatmap_interferon.png", bg = "white")


##
### 180mV vs 100mV
##
# perform an exact test for 180mV vs 100mV
tested <- exactTest(list, pair=c("100mV", "180mV"))

# create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table %>% rownames_to_column(var="gene")


# write results table of DE genes to a CSV file
write.csv(resultsTbl, file="180mV_100mV_topTags.csv")


# plot log-fold change against log-counts per million, with DE genes highlighted
# the blue lines indicate 2-fold changes
jpeg("180mV_100mV_plotMD.jpg")
plotMD(tested)
abline(h=c(-1, 1), col="grey")
dev.off()


# identify significantly DE genes
resultsTbl$sigDE <- "NA"
resultsTbl$sigDE[resultsTbl$logFC > 1 & resultsTbl$FDR < 0.05] <- "UP"
resultsTbl$sigDE[resultsTbl$logFC < -1 & resultsTbl$FDR < 0.05] <- "DOWN"


# create volcano plot
jpeg("180mV_100mV_plotVolcano.jpg")
ggplot(data=resultsTbl, aes(x=logFC, y=-log10(FDR), color = sigDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = color_subset)
dev.off()


# remove the sigDE column
resultsTbl_clean <- resultsTbl[ , -which(names(resultsTbl) %in% c("sigDE"))]


# return all rows from resultsTbl_clean and genes_att
resultsTbl_att <- resultsTbl_clean %>% left_join(genes_att, by = c("gene" = "ensembl_gene_id"))


# sort the results by FDR
resultsTbl_ann <- resultsTbl_att[order(resultsTbl_att$FDR),]


# write the sorted DE genes to a CSV file
write.csv(resultsTbl_ann, file="180mV_100mV_annotated.csv")


