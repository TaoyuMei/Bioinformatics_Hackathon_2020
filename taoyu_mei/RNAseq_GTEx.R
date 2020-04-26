{
library(R.utils)
library(tidyverse)
library(DESeq2)
library(readxl)
# library(cmapR)
library(stringr)
library(ggplot2)
library(clusterProfiler)
}

# setwd("D:/GitHub/Bioinformatics_Hackathon_2020/taoyu_mei")


# download preprocessed RNA-Seq data from GTEx ----------------------------

# download.file("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz", 
#               destfile = "./GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")
# R.utils::gunzip("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")

geneCounts <- read_tsv("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct",
                       skip = 2)

# geneCounts <- read.gct("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
# geneCounts <- cmapR::parse_gctx("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")


# subset the count matrix by patients included ---------------------------

# pheno1 <- read_excel("../simons_eda/GTEx_pancreas_liver_images_liverfat_pancreasfat.xlsx")
pheno2 <- read_excel("../simons_eda/GTEx_pancreas_liver_images_liverfat_pancreasfat_seq.xlsx")
# the first file include all patients with images thus fat percentage
# the second file only include patients with RNA-Seq data, so it's a subset of the first

# pheno1$Tissue.Sample.ID_pancreas[1] %in% colnames(geneCounts)

# SAMID <- c(pheno2$SAMPID_pancreas, pheno2$SAMPID_liver)
# SAMID <- intersect(SAMID, colnames(geneCounts))
# geneCounts_subsetted <- geneCounts[c('Name', 'Description', SAMID)]


# download.file("https://www.dropbox.com/s/0ezaxy5e56nbp34/count_matrix_target_subset.tsv.gz",
#               destfile = "count_matrix_target_subset.tsv.gz")
# R.utils::gunzip("count_matrix_target_subset.tsv.gz")
# geneCounts_subsetted <- read_tsv("count_matrix_target_subset.tsv")

geneCounts_pancreas <- geneCounts[c('Name', 'Description', 
                                    intersect(pheno2$SAMPID_pancreas, colnames(geneCounts)))]

geneCounts_liver <- geneCounts[c('Name', 'Description', 
                                 intersect(pheno2$SAMPID_liver, colnames(geneCounts)))]



# make ready for differential expression analysis -------------------------

# make 2 colData for the 2 organs; add a group variable according to age group

# pancreas
colDataPancreas <- pheno2[c('SAMPID_pancreas', 'Sex', 'Age.Bracket', 'Hardy.Scale', 
                            'Fat,Percentage_liver', 'Fat,Percentage_pancreas')]
colDataPancreas <- filter(colDataPancreas, 
                          SAMPID_pancreas %in% colnames(geneCounts_pancreas))
colDataPancreas <- unique.data.frame(colDataPancreas)  # remove duplicate

colDataPancreas$Age.Bracket <- gsub("^[60|70].*", "old", colDataPancreas$Age.Bracket)
colDataPancreas$Age.Bracket <- gsub("^[40|50].*", "mid_age", colDataPancreas$Age.Bracket)
colDataPancreas$Age.Bracket <- gsub("^[20|30].*", "young", colDataPancreas$Age.Bracket)

row_names <- colDataPancreas$SAMPID_pancreas
colDataPancreas <-as.matrix(colDataPancreas[-1])
row.names(colDataPancreas) <- row_names


# liver
colDataLiver <- pheno2[c('SAMPID_liver', 'Sex', 'Age.Bracket', 'Hardy.Scale', 
                            'Fat,Percentage_liver', 'Fat,Percentage_pancreas')]
colDataLiver <- filter(colDataLiver, 
                          SAMPID_liver %in% colnames(geneCounts_liver))
colDataLiver <- unique.data.frame(colDataLiver)

colDataLiver$Age.Bracket <- gsub("^[60|70].*", "old", colDataLiver$Age.Bracket)
colDataLiver$Age.Bracket <- gsub("^[40|50].*", "mid_age", colDataLiver$Age.Bracket)
colDataLiver$Age.Bracket <- gsub("^[20|30].*", "young", colDataLiver$Age.Bracket)

row_names <- colDataLiver$SAMPID_liver
colDataLiver <-as.matrix(colDataLiver[-1])
row.names(colDataLiver) <- row_names


# convert the count matrices' row names to fit DESeq2
countMatrixPancreas <- as.matrix(geneCounts_pancreas[c(-1, -2)])
row.names(countMatrixPancreas) <- geneCounts_pancreas$Name

countMatrixLiver <- as.matrix(geneCounts_liver[c(-1, -2)])
row.names(countMatrixLiver) <- geneCounts_liver$Name


# check the sample existence and order

# all(rownames(colDataLiver) %in% colnames(countMatrixLiver))
all(rownames(colDataLiver) == colnames(countMatrixLiver))
all(rownames(colDataPancreas) == colnames(countMatrixPancreas))


# re-order according to the variable of interest (control level at first)
colDataPancreas <-colDataPancreas[order(colDataPancreas[, "Age.Bracket"], 
                                        decreasing = TRUE), ]
countMatrixPancreas <- countMatrixPancreas[, rownames(colDataPancreas)]

colDataLiver <-colDataLiver[order(colDataLiver[, "Age.Bracket"], 
                                  decreasing = TRUE), ]
countMatrixLiver <- countMatrixLiver[, rownames(colDataLiver)]


# output the matrices for other pipelines
saveRDS(colDataLiver, "colDataLiver.rds")
saveRDS(countMatrixLiver, "countMatrixLiver.rds")
saveRDS(colDataPancreas, "colDataPancreas.rds")
saveRDS(countMatrixPancreas, "countMatrixPancreas.rds")

write_tsv(as.data.frame(colDataLiver), "colDataLiver.tsv")
write_tsv(as.data.frame(countMatrixLiver), "countMatrixLiver.tsv")
write_tsv(as.data.frame(colDataPancreas), "colDataPancreas.tsv")
write_tsv(as.data.frame(countMatrixPancreas), "countMatrixPancreas.tsv")

# load from RDS on the next day to save time

{
colDataLiver <- readRDS("colDataLiver.rds")
countMatrixLiver <- readRDS("countMatrixLiver.rds")
colDataPancreas <- readRDS("colDataPancreas.rds")
countMatrixPancreas <- readRDS("countMatrixPancreas.rds")
}

# create a dds object
ddsPancreas <- DESeqDataSetFromMatrix(countData = countMatrixPancreas, 
                                      colData = colDataPancreas,
                                      design = ~ Sex + Hardy.Scale + Age.Bracket)


ddsLiver <- DESeqDataSetFromMatrix(countData = countMatrixLiver, 
                                      colData = colDataLiver,
                                      design = ~ Sex + Hardy.Scale + Age.Bracket)
# column names with ',' cause errors



# overview of the gene expression before DE analysis ----------------------
vsdPancreas <- vst(ddsPancreas)
vsdLiver <- vst(ddsLiver)

# PCA
pcaDataPancreas <- DESeq2::plotPCA(vsdPancreas, intgroup = c("Sex", "Age.Bracket"),
                                   returnData = TRUE)
ggplot(data = pcaDataPancreas, aes(PC1, PC2, color = Sex, shape = Age.Bracket)) +
  geom_point(size=3)

pcaDataLiver <- DESeq2::plotPCA(vsdLiver, intgroup = c("Sex", "Age.Bracket"),
                                returnData = TRUE)
ggplot(pcaDataLiver, aes(PC1, PC2, color = Sex, shape = Age.Bracket)) +
  geom_point(size=3)

# obviously, gender contribute much more than age to the variance in gene expression
# gender need to be adjusted when performing DE analysis between age groups.


# heatmap 


# hierarchical clustering tree 



# DE analysis -------------------------------------------------------------

ddsPancreas <- DESeq(ddsPancreas)
resultsNames(ddsPancreas) # lists the coefficients 
###### why only "Age.Bracket_old_vs_mid_age" and "Age.Bracket_young_vs_mid_age"
###### but no "Age.Bracket_young_vs_old" ?

res <- results(ddsPancreas, name="")

ddsLiver <- DESeq()




# clusterProfiler (GSEA etc.) ---------------------------------------------



