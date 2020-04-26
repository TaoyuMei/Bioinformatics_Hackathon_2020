#Bioneer is describing combined effect of variation i.e. clonal variation and gwoth phase specific variation
#here we have followed the same description as metabolic growth phase only by excluding the parantal cell lines fed with hGL
BiocManager::install("CePa")
BiocManager::install("Piano")
BiocManager::install("xlsx")
library(CePa)
require(DESeq2)
library(edgeR)
library(limma)
#library(Glimma) #package ‘Glimma’ is not available (for R version 3.3.0)
library(gplots)
library(dplyr)
#library(org.Mm.eg.db)
#library(org.Hs.eg.db)
library(RColorBrewer)
#Read all of the htseq files woth the help of as.matrix

head(keep)
counts.keep <- col[keep,] # Subset the rows of countdata to keep the more highly expressed genes
dim(counts.keep)  
summary(keep) # This command will summerize the gene sets passed filter of 0.5. And total of 16032 genes were having CPM count more than .5 as it was found to be true and 10603 genes failed to pass this filter.  

colDataLiver <- readRDS(file = "colDataLiver.rds")
colDataPancreas <- readRDS(file = "colDataPancreas.rds")
countMatrixLiver <- readRDS(file = "countMatrixLiver.rds")
countMatrixPancreas <- readRDS(file = "countMatrixPancreas.rds")

colnames(countMatrixLiver)
y <- DGEList(count = countMatrixLiver, group = factor)

#y <- DGEList(count = counts.keep, group = group) # DGE list object in edgeR is used to store count data
ynorm <- calcNormFactors(y)

factor
design <- model.matrix(~ 0 + factor)
design
colnames(design) <- levels(factor)
design


v <- voom(ynorm,design,plot = TRUE) #This plot can also tell us if there are any genes that look really variable in our data, and if we’ve filtered the low counts adequately.
#voom: transforms the readcounts into logCPM while taking into consideration the mean variance of the data
v
names(v)
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(v$E),col="blue")
fit <- lmFit(v,design) # lmFit function was used to find out the differentially expresssed genes. It uses voom transformed function along with the design matix that we have already specified in the command.
#fit <- glmQLFit(v)
#lmFit calculates group's mean based on the design matirix as well as gene wise variance.
names(fit)

#this contrast is basically on the basis of growth phase varoation 
cont.matrix1 <- makeContrasts(Clone1 = young-old,levels=design)
cont.matrix2 <- makeContrasts(Clone2 = mid_age-old,levels=design)
cont.matrix3 <- makeContrasts(Clone3 = young-mid_age,levels=design)
#NR is indicating noise reduction 
#This is only for vertical contrast 


fit.cont1 <- contrasts.fit(fit, cont.matrix1) #Then we fit the contrast matrix with the help of fit function 
fit.cont1 <- eBayes(fit.cont1)
dim(fit.cont1)
summa.fit1 <- decideTests(fit.cont1)
summary(summa.fit1)
Contrast_Clone1 <- topTable(fit.cont1,coef="Clone1",sort.by="p", number = "Inf")
Contrast_Clone1 <- cbind(rownames(Contrast_Clone1), Contrast_Clone1)
rownames(Contrast_Clone1) <- NULL
colnames(Contrast_Clone1) <- c("ENSEMBL", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")

Contrast_Clone1$ENSEMBL=gsub("\\..*","",Contrast_Clone1$ENSEMBL)
gene.df1 <- bitr(Contrast_Clone1$ENSEMBL, fromType = "ENSEMBL", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db)
Contrast_Clone1_new <- merge(gene.df1, Contrast_Clone1, by=c("ENSEMBL"))
write.table(Contrast_Clone1_new,"Clone1_new.txt")


fit.cont2 <- contrasts.fit(fit, cont.matrix2) #Then we fit the contrast matrix with the help of fit function 
fit.cont2 <- eBayes(fit.cont2)
dim(fit.cont2)
summa.fit2 <- decideTests(fit.cont2)
summary(summa.fit2)
Contrast_Clone2 <- topTable(fit.cont2,coef="Clone2",sort.by="p", number = "Inf")
Contrast_Clone2 <- cbind(rownames(Contrast_Clone2), Contrast_Clone2)
rownames(Contrast_Clone2) <- NULL
colnames(Contrast_Clone2) <- c("ENSEMBL", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")

Contrast_Clone2$ENSEMBL=gsub("\\..*","",Contrast_Clone2$ENSEMBL)
gene.df2 <- bitr(Contrast_Clone2$ENSEMBL, fromType = "ENSEMBL", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db)
Contrast_Clone2_new <- merge(gene.df2, Contrast_Clone2, by=c("ENSEMBL"))
write.table(Contrast_Clone2_new,"Clone2_new.txt")

fit.cont3 <- contrasts.fit(fit, cont.matrix3) #Then we fit the contrast matrix with the help of fit function 
fit.cont3 <- eBayes(fit.cont3)
dim(fit.cont3)
summa.fit3 <- decideTests(fit.cont3)
summary(summa.fit3)
Contrast_Clone3 <- topTable(fit.cont3,coef="Clone3",sort.by="p", number = "Inf")
#row.names(Contrast_Clone3)=gsub("\\..*","",row.names(Contrast_Clone3))
Contrast_Clone3 <- cbind(rownames(Contrast_Clone3), Contrast_Clone3)
rownames(Contrast_Clone3) <- NULL
colnames(Contrast_Clone3) <- c("ENSEMBL", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")

Contrast_Clone3$ENSEMBL=gsub("\\..*","",Contrast_Clone3$ENSEMBL)
gene.df3 <- bitr(Contrast_Clone3$ENSEMBL, fromType = "ENSEMBL", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db)
Contrast_Clone3_new <- merge(gene.df3, Contrast_Clone3, by=c("ENSEMBL"))
write.table(Contrast_Clone3_new,"Clone3_new.txt")
