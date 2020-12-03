## prepare datasets and conduct eQTL analysis
## batch effect adjustment for gene expression
## filter SNPs: MAF>0.01, Info score>=0.7
## D Huo 

 rm(list = ls())
 #ls()
 con <- file("/gpfs/data/huo-lab/ROOT/Tunde/eQTL/codes/eQTL_analysis2.R.log")
 sink(con) 

 setwd("/gpfs/data/huo-lab/ROOT/Tunde/eQTL") 
 ## setwd("Y:/ROOT/Tunde/eQTL") 
 
#### Section [1] genotype files: geno, covariates
## sample file
sample.geno <- read.table("TCGA_MutSig.sample_v1.sample", sep = " ", header = TRUE, stringsAsFactors = FALSE) 

dim(sample.geno)
head(sample.geno[, c("id_1", "age", "er", "pr", "her2", "race_geno", "race_repo", "ev1")])
head(sample.geno[-1, c("id_1", "age", "er", "pr", "her2", "race_geno", "race_repo", "ev1")])
sample.geno1 <- sample.geno[-1, c("id_1", "age", "er", "pr", "her2", "race_geno", "race_repo", "prop_ceu", "prop_yri", "prop_asn", "ev1", "ev2", "ev3", "ev4", "ev5", "ev6", "ev7", "ev8", "ev9", "ev10")]
dim(sample.geno1)
table(sample.geno1$race_geno)
sample.geno1$race_geno[is.na(sample.geno1$race_geno)] <- 99 
table(sample.geno1$race_geno)

list.B <- sample.geno1$race_geno==1 
list.W <- sample.geno1$race_geno==2
table(list.B)
table(list.W)
sample.geno.B <- sample.geno1[list.B, ]
sample.geno.W <- sample.geno1[list.W, ]

## imputation info file of 4 chr:  
setwd("/gpfs/data/huo-lab/TCGA/SNP6_level2/Gao/imputed_data_concatenate_by_chrom_info") 
## setwd("Y:/TCGA/SNP6_level2/Gao/imputed_data_concatenate_by_chrom_info") 
info1 <- read.table("tcga.chr1.phased.impute2.info", sep = " ", header = TRUE, stringsAsFactors = FALSE, as.is=TRUE)
info1A <- info1[(info1$position>=109679756 & info1$position<=110679756) | (info1$position>=214830292 & info1$position<=215830292), c(2:5,7)]

info5 <- read.table("tcga.chr5.phased.impute2.info", sep = " ", header = TRUE, stringsAsFactors = FALSE, as.is=TRUE)
info5A <- info5[(info5$position>=131301947  & info5$position<=132301947 ), c(2:5,7)]

info7 <- read.table("tcga.chr7.phased.impute2.info", sep = " ", header = TRUE, stringsAsFactors = FALSE, as.is=TRUE)
info7A <- info7[(info7$position>=73859358   & info7$position<=74859358  ), c(2:5,7)]

info15 <- read.table("tcga.chr15.phased.impute2.info", sep = " ", header = F, stringsAsFactors = FALSE, as.is=TRUE)
colnames(info15) <- colnames(info1)
info15A <- info15[(info15$position>=74647332 & info15$position<=76228474) | (info15$position>=100407094 & info15$position<=101407094), c(2:5,7)]

rm("info1", "info5", "info7", "info15")


## genotypes of 6 regions
setwd("/gpfs/data/huo-lab/ROOT/Tunde/eQTL") 
## setwd("Y:/ROOT/Tunde/eQTL") 

geno.1p13 <- read.table("TCGA_1p13.3.gen.gz", sep = " ", header = FALSE, stringsAsFactors = FALSE, as.is=TRUE) 
dim(geno.1p13)
geno.1p13[1:3, 1:10]

geno.1p13A <-  merge(geno.1p13, info1A, by.x="V2", by.y="rs_id", all.x=T)
dim(geno.1p13A)
summary(geno.1p13A[, c("exp_freq_a1", "info", "type")]) 
summary(geno.1p13A$exp_freq_a1<0.95 & geno.1p13A$exp_freq_a1>0.05)
summary(geno.1p13A$exp_freq_a1<0.99 & geno.1p13A$exp_freq_a1>0.01) 
geno.1p13B <- geno.1p13A[geno.1p13A$exp_freq_a1<0.99 & geno.1p13A$exp_freq_a1>0.01, ]

geno.1p13.left <- geno.1p13B[, c(1:5,3170:3172)]   
geno.1p13.matrix <- as.matrix(geno.1p13B[, 6:3167]) 

dose.1p13 <- matrix(data=NA, nrow=nrow(geno.1p13.left), ncol=1054) 
for (j in (1:1054)) {  
   d = geno.1p13.matrix[,j*3-1] + geno.1p13.matrix[,j*3]*2
   dose.1p13[,j]<-d
}

dose.1p13.B <- dose.1p13[, list.B]
dose.1p13.W <- dose.1p13[, list.W]
rownames(dose.1p13.B) <- geno.1p13.left$V2
rownames(dose.1p13.W) <- geno.1p13.left$V2
dose.1p13.B[1:3, 1:3]
dose.1p13.W[1:3, 1:3]


geno.1q41 <- read.table("TCGA_1q41.gen.gz", sep = " ", header = FALSE, stringsAsFactors = FALSE, as.is=TRUE) 
dim(geno.1q41)
geno.1q41[1:3, 1:10]

geno.1q41A <-  merge(geno.1q41, info1A, by.x="V2", by.y="rs_id", all.x=T)
dim(geno.1q41A)
summary(geno.1q41A[, c("exp_freq_a1", "info", "type")]) 
summary(geno.1q41A$exp_freq_a1<0.95 & geno.1q41A$exp_freq_a1>0.05)
summary(geno.1q41A$exp_freq_a1<0.99 & geno.1q41A$exp_freq_a1>0.01) 
geno.1q41B <- geno.1q41A[geno.1q41A$exp_freq_a1<0.99 & geno.1q41A$exp_freq_a1>0.01, ]

geno.1q41.left <- geno.1q41B[, c(1:5,3170:3172)]
geno.1q41.matrix <- as.matrix(geno.1q41B[, 6:3167]) 

dose.1q41 <- matrix(, nrow=nrow(geno.1q41.left), ncol=1054) 
for (j in (1:1054)) {  
   d = geno.1q41.matrix[,j*3-1] + geno.1q41.matrix[,j*3]*2
   dose.1q41[,j]<-d
}

dose.1q41.B <- dose.1q41[, list.B]
dose.1q41.W <- dose.1q41[, list.W]
rownames(dose.1q41.B) <- geno.1q41.left$V2
rownames(dose.1q41.W) <- geno.1q41.left$V2
dose.1q41.B[1:3, 1:3]
dose.1q41.W[1:3, 1:3]


geno.5q31 <- read.table("TCGA_5q31.1.gen.gz", sep = " ", header = FALSE, stringsAsFactors = FALSE, as.is=TRUE) 
dim(geno.5q31)
geno.5q31[1:3, 1:10]

geno.5q31A <-  merge(geno.5q31, info5A, by.x="V2", by.y="rs_id", all.x=T)
dim(geno.5q31A)
summary(geno.5q31A[, c("exp_freq_a1", "info", "type")]) 
summary(geno.5q31A$exp_freq_a1<0.95 & geno.5q31A$exp_freq_a1>0.05)
summary(geno.5q31A$exp_freq_a1<0.99 & geno.5q31A$exp_freq_a1>0.01) 
geno.5q31B <- geno.5q31A[geno.5q31A$exp_freq_a1<0.99 & geno.5q31A$exp_freq_a1>0.01, ]

geno.5q31.left <- geno.5q31B[, c(1:5,3170:3172)]
geno.5q31.matrix <- as.matrix(geno.5q31B[, 6:3167]) 

dose.5q31 <- matrix(, nrow=nrow(geno.5q31.left), ncol=1054) 
for (j in (1:1054)) {  
   d = geno.5q31.matrix[,j*3-1] + geno.5q31.matrix[,j*3]*2
   dose.5q31[,j]<-d
}

dose.5q31.B <- dose.5q31[, list.B]
dose.5q31.W <- dose.5q31[, list.W]
rownames(dose.5q31.B) <- geno.5q31.left$V2
rownames(dose.5q31.W) <- geno.5q31.left$V2
dose.5q31.B[1:3, 1:3]
dose.5q31.W[1:3, 1:3]


geno.7q11 <- read.table("TCGA_7q11.23.gen.gz", sep = " ", header = FALSE, stringsAsFactors = FALSE, as.is=TRUE) 
dim(geno.7q11)
geno.7q11[1:3, 1:10]

geno.7q11A <-  merge(geno.7q11, info7A, by.x="V2", by.y="rs_id", all.x=T)
dim(geno.7q11A)
summary(geno.7q11A[, c("exp_freq_a1", "info", "type")]) 
summary(geno.7q11A$exp_freq_a1<0.95 & geno.7q11A$exp_freq_a1>0.05)
summary(geno.7q11A$exp_freq_a1<0.99 & geno.7q11A$exp_freq_a1>0.01) 
geno.7q11B <- geno.7q11A[geno.7q11A$exp_freq_a1<0.99 & geno.7q11A$exp_freq_a1>0.01, ]

geno.7q11.left <- geno.7q11B[, c(1:5,3170:3172)]
geno.7q11.matrix <- as.matrix(geno.7q11B[, 6:3167]) 

dose.7q11 <- matrix(data=NA, nrow=nrow(geno.7q11.left), ncol=1054) 
for (j in (1:1054)) {  
   d = geno.7q11.matrix[,j*3-1] + geno.7q11.matrix[,j*3]*2
   dose.7q11[,j]<-d
}

dose.7q11.B <- dose.7q11[, list.B]
dose.7q11.W <- dose.7q11[, list.W]
rownames(dose.7q11.B) <- geno.7q11.left$V2
rownames(dose.7q11.W) <- geno.7q11.left$V2
dose.7q11.B[1:3, 1:3]
dose.7q11.W[1:3, 1:3]
   

geno.15q24 <- read.table("TCGA_15q24.gen.gz", sep = " ", header = FALSE, stringsAsFactors = FALSE, as.is=TRUE) 
dim(geno.15q24)
geno.15q24[1:3, 1:10]

geno.15q24A <-  merge(geno.15q24, info15A, by.x="V2", by.y="rs_id", all.x=T)
dim(geno.15q24A)
summary(geno.15q24A[, c("exp_freq_a1", "info", "type")]) 
summary(geno.15q24A$exp_freq_a1<0.95 & geno.15q24A$exp_freq_a1>0.05)
summary(geno.15q24A$exp_freq_a1<0.99 & geno.15q24A$exp_freq_a1>0.01) 
geno.15q24B <- geno.15q24A[geno.15q24A$exp_freq_a1<0.99 & geno.15q24A$exp_freq_a1>0.01, ]

geno.15q24.left <- geno.15q24B[, c(1:5,3170:3172)]
geno.15q24.matrix <- as.matrix(geno.15q24B[, 6:3167]) 

dose.15q24 <- matrix(data=NA, nrow=nrow(geno.15q24.left), ncol=1054) 
for (j in (1:1054)) {  
   d = geno.15q24.matrix[,j*3-1] + geno.15q24.matrix[,j*3]*2
   dose.15q24[,j]<-d
}

dose.15q24.B <- dose.15q24[, list.B]
dose.15q24.W <- dose.15q24[, list.W]
rownames(dose.15q24.B) <- geno.15q24.left$V2
rownames(dose.15q24.W) <- geno.15q24.left$V2
dose.15q24.B[1:3, 1:3]
dose.15q24.W[1:3, 1:3]


geno.15q26 <- read.table("TCGA_15q26.3.gen.gz", sep = " ", header = FALSE, stringsAsFactors = FALSE, as.is=TRUE) 
dim(geno.15q26)
geno.15q26[1:3, 1:10]

geno.15q26A <-  merge(geno.15q26, info15A, by.x="V2", by.y="rs_id", all.x=T)
dim(geno.15q26A)
summary(geno.15q26A[, c("exp_freq_a1", "info", "type")]) 
summary(geno.15q26A$exp_freq_a1<0.95 & geno.15q26A$exp_freq_a1>0.05)
summary(geno.15q26A$exp_freq_a1<0.99 & geno.15q26A$exp_freq_a1>0.01) 
geno.15q26B <- geno.15q26A[geno.15q26A$exp_freq_a1<0.99 & geno.15q26A$exp_freq_a1>0.01, ]

geno.15q26.left <- geno.15q26B[, c(1:5,3170:3172)]
geno.15q26.matrix <- as.matrix(geno.15q26B[, 6:3167]) 

dose.15q26 <- matrix(data=NA, nrow=nrow(geno.15q26.left), ncol=1054) 
for (j in (1:1054)) {  
   d = geno.15q26.matrix[,j*3-1] + geno.15q26.matrix[,j*3]*2
   dose.15q26[,j]<-d
}

dose.15q26.B <- dose.15q26[, list.B]
dose.15q26.W <- dose.15q26[, list.W]
rownames(dose.15q26.B) <- geno.15q26.left$V2
rownames(dose.15q26.W) <- geno.15q26.left$V2
dose.15q26.B[1:3, 1:3]
dose.15q26.W[1:3, 1:3]




#### Section [2] gene expression files and gene list
## gene expression
setwd("/gpfs/data/huo-lab/TCGA/RNASeqV2/Matrix") 
## setwd("Y:/TCGA/RNASeqV2/Matrix") 

## gene.express <- read.table("matrix_log2_18811genes_928samples.txt", sep="\t", row.names=1, header = TRUE, stringsAsFactors = FALSE) 
load("/gpfs/data/huo-lab/TCGA/RNASeqV2/Matrix/RNAseq_log2_18811genes_928samples.RData")
ls()
dim(BWlog)
BWlog[1:3, 1:6]
dim(pheno1)
head(pheno1)

## most but not all gene symbols are in upper case. 
gene.name <- data.frame(matrix(unlist(strsplit(rownames(BWlog), '[|]')), ncol = 2, byrow = TRUE), stringsAsFactors = FALSE)
colnames(gene.name) <- c("genesymbol", "GeneID")
gene.name$order <- 1:nrow(gene.name)


setwd("/gpfs/data/huo-lab/ROOT/Tunde/eQTL") 
## setwd("Y:/ROOT/Tunde/eQTL") 
gene.list <- read.table("Gene_position_6Loci.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) 

gene.comb <- merge(gene.list, gene.name, by.x="gene_name", by.y="genesymbol", all=TRUE)
head(gene.comb)
table(is.na(gene.comb$order), is.na(gene.comb$gene_biotype))
## not merged: 
table(gene.comb$gene_biotype[is.na(gene.comb$order)]) 
## list of not merged protein coding genes
gene.comb[is.na(gene.comb$order) & gene.comb$gene_biotype=="protein_coding" , ]
## merged: 
table(gene.comb$gene_biotype[!is.na(gene.comb$order)]) 

gene.comb1 <- gene.comb[!is.na(gene.comb$order),]
gene.comb2 <- gene.comb1[order(gene.comb1$order),]
gene.anno <- gene.comb2[!is.na(gene.comb2$loci), ]
   
GeneExp <- BWlog[!is.na(gene.comb2$loci), ]
dim(GeneExp) 
GeneExp[1:5, 1:5]

table(pheno1$RACEgeno)
pheno1.B <- pheno1[pheno1$RACEgeno=="Black", ]
GeneExp.B <- GeneExp[, pheno1$RACEgeno=="Black"]

pheno1.W <- pheno1[pheno1$RACEgeno=="White", ]
GeneExp.W <- GeneExp[, pheno1$RACEgeno=="White"]


#### Section [3] CNV data
setwd("/gpfs/data/huo-lab/ROOT/Tunde/eQTL/CNV_eqtl") 
## setwd("Y:/ROOT/Tunde/eQTL/CNV_eqtl") 

CNV.B.H <-  read.table("all_thresholded.by_genes_Afr_Amer.txt", sep = "\t", header = FALSE, as.is=TRUE, stringsAsFactors = FALSE) 
CNV.W.H <-  read.table("all_thresholded.by_genes_caucasian.txt", sep = "\t", header = FALSE, as.is=TRUE, stringsAsFactors = FALSE) 
CNV.B.H1 <- t(as.matrix(CNV.B.H[1, ])) 
CNV.B.H2 <- substr(CNV.B.H1, 1, 12)
CNV.B.H3 <- sub(" ", "", CNV.B.H2)
CNV.W.H1 <- t(as.matrix(CNV.W.H[1, ])) 
CNV.W.H2 <- substr(CNV.W.H1, 1, 12)
CNV.W.H3 <- sub(" ", "", CNV.W.H2)
rm("CNV.B.H", "CNV.W.H")

CNV.B <-  read.table("all_thresholded.by_genes_Afr_Amer.txt", sep = "\t", header = TRUE, as.is=TRUE, stringsAsFactors = FALSE) 
CNV.W <-  read.table("all_thresholded.by_genes_caucasian.txt", sep = "\t", header = TRUE, as.is=TRUE, stringsAsFactors = FALSE) 
CNV.B[1:5, 1:5]
colnames(CNV.B) <- CNV.B.H3
colnames(CNV.W) <- CNV.W.H3

CNV.B.data <- CNV.B[CNV.B$GeneSymbol %in% gene.anno$gene_name, ] 
CNV.B.data2 <-CNV.B.data[order(CNV.B.data$GeneSymbol),] 
CNV.B.data3 <- matrix(data=NA, nrow=nrow(CNV.B.data2), ncol=(ncol(CNV.B.data2)-3)) 
d <- vector(mode = "numeric", length=nrow(CNV.B.data2))
for (j in (4:ncol(CNV.B.data2))) {  
   d[CNV.B.data2[, j]== -2] <- 0 
   d[CNV.B.data2[, j]> -2 & CNV.B.data2[, j]< 2] <- 1 
   d[CNV.B.data2[, j]== 2] <- 2
   CNV.B.data3[, j-3] <- t(d)
}
CNV.B.data4 <- cbind(CNV.B.data2[, 1:3], CNV.B.data3) 
colnames(CNV.B.data4) = colnames(CNV.B.data2)

CNV.W.data <- CNV.W[CNV.W$GeneSymbol %in% gene.anno$gene_name, ] 
CNV.W.data2 <-CNV.W.data[order(CNV.W.data$GeneSymbol),] 
CNV.W.data3 <- matrix(data=NA, nrow=nrow(CNV.W.data2), ncol=(ncol(CNV.W.data2)-3))  
d <- vector(mode = "numeric", length=nrow(CNV.W.data2))
for (j in (4:ncol(CNV.W.data2))) {  
   d[CNV.W.data2[, j]== -2] <- 0 
   d[CNV.W.data2[, j]> -2 & CNV.W.data2[, j]< 2] <- 1 
   d[CNV.W.data2[, j]== 2] <- 2
   CNV.W.data3[, j-3] <- t(d)
}
CNV.W.data4 <- cbind(CNV.W.data2[, 1:3], CNV.W.data3) 
colnames(CNV.W.data4) = colnames(CNV.W.data2)
dim(CNV.W.data4)


## section [4] merge data together

dim(sample.geno.B)
sample.geno.B$col1 <- 1:nrow(sample.geno.B)
   
ID.B <- colnames(CNV.B.data4)[4:ncol(CNV.B.data4)]
order3 <- 1:length(ID.B)
sample.CNV.B <- cbind.data.frame(ID.B, order3) 

dim(pheno1.B) 
pheno1.B$col2 <- 1:nrow(pheno1.B)

sample.B.1 <- merge(sample.geno.B, pheno1.B, by.x="id_1", by.y="patientid", all=TRUE) 
table(is.na(sample.B.1$col1), is.na(sample.B.1$col2))

sample.B.2 <- merge(sample.B.1, sample.CNV.B, by.x="id_1", by.y="ID.B", all=TRUE) 
table(is.na(sample.B.2$col1), is.na(sample.B.2$col2), is.na(sample.B.2$order3))

sample.B.3 <- sample.B.2[!is.na(sample.B.2$col1) & !is.na(sample.B.2$col2) & !is.na(sample.B.2$order3), ]
sample.B.3 <- sample.B.3[order(sample.B.3$id_1), ] 
sample.B.3$orderF <- 1:nrow(sample.B.3)

## merge back to genotype file
sample.geno.B1 <- merge(sample.geno.B, sample.B.3, by.x="id_1", by.y="id_1", all=TRUE) 
sample.geno.B1 <- sample.geno.B1[order(sample.geno.B1$col1.x) , ] 

sample.geno.B2 <- sample.geno.B1[!is.na(sample.geno.B1$orderF) , ]

   b <- dose.7q11.B[, !is.na(sample.geno.B1$orderF)] 
   dose.7q11.B.F <- b[, order(sample.geno.B2$orderF)]
   
   b <- dose.5q31.B[, !is.na(sample.geno.B1$orderF)] 
   dose.5q31.B.F <- b[, order(sample.geno.B2$orderF)]
 
   b <- dose.15q26.B[, !is.na(sample.geno.B1$orderF)] 
   dose.15q26.B.F <- b[, order(sample.geno.B2$orderF)]
   
   b <- dose.1p13.B[, !is.na(sample.geno.B1$orderF)] 
   dose.1p13.B.F <- b[, order(sample.geno.B2$orderF)]

   b <- dose.15q24.B[, !is.na(sample.geno.B1$orderF)] 
   dose.15q24.B.F <- b[, order(sample.geno.B2$orderF)]

   b <- dose.1q41.B[, !is.na(sample.geno.B1$orderF)] 
   dose.1q41.B.F <- b[, order(sample.geno.B2$orderF)]
   
## merge back to gen expression file   
   pheno1.B1 <- merge(pheno1.B, sample.B.3, by.x="patientid", by.y="id_1", all=TRUE) 
   pheno1.B1 <- pheno1.B1[order(pheno1.B1$col2.x) , ] 
   pheno1.B2 <- pheno1.B1[!is.na(pheno1.B1$orderF) , ]
   
   tmp <- GeneExp.B[, !is.na(pheno1.B1$orderF)] 
   GeneExp.B.F <- tmp[, order(pheno1.B2$orderF)]
   
## merge back to CNV file      
   sample.CNV.B1 <- merge(sample.CNV.B, sample.B.3, by.x="ID.B", by.y="id_1", all=TRUE) 
   sample.CNV.B1 <- sample.CNV.B1[order(sample.CNV.B1$order3.x) , ] 
   sample.CNV.B2 <- sample.CNV.B1[!is.na(sample.CNV.B1$orderF) , ]
   
   CNV.B2 <- CNV.B.data4[, 4:ncol(CNV.B.data4)]
   
   tmp <- CNV.B2[, !is.na(sample.CNV.B1$orderF)] 
   CNV.B.F <- tmp[, order(sample.CNV.B2$orderF)]
   rownames(CNV.B.F) <- CNV.B.data4$GeneSymbol
   
   
## for whites
   dim(sample.geno.W)
   sample.geno.W$col1 <- 1:nrow(sample.geno.W)
   
   ID.W <- colnames(CNV.W.data4)[4:ncol(CNV.W.data4)]
   order3 <- 1:length(ID.W)
   sample.CNV.W <- cbind.data.frame(ID.W, order3) 
   
   dim(pheno1.W) 
   pheno1.W$col2 <- 1:nrow(pheno1.W)
   
   sample.W.1 <- merge(sample.geno.W, pheno1.W, by.x="id_1", by.y="patientid", all=TRUE) 
   table(is.na(sample.W.1$col1), is.na(sample.W.1$col2))
   
   sample.W.2 <- merge(sample.W.1, sample.CNV.W, by.x="id_1", by.y="ID.W", all=TRUE) 
   table(is.na(sample.W.2$col1), is.na(sample.W.2$col2), is.na(sample.W.2$order3))
   
   sample.W.3 <- sample.W.2[!is.na(sample.W.2$col1) & !is.na(sample.W.2$col2) & !is.na(sample.W.2$order3), ]
   sample.W.3 <- sample.W.3[order(sample.W.3$id_1), ] 
   sample.W.3$orderF <- 1:nrow(sample.W.3)
   
   ## merge back to genotype file
   sample.geno.W1 <- merge(sample.geno.W, sample.W.3, by.x="id_1", by.y="id_1", all=TRUE) 
   sample.geno.W1 <- sample.geno.W1[order(sample.geno.W1$col1.x) , ] 
   
   sample.geno.W2 <- sample.geno.W1[!is.na(sample.geno.W1$orderF) , ]
   
   b <- dose.7q11.W[, !is.na(sample.geno.W1$orderF)] 
   dose.7q11.W.F <- b[, order(sample.geno.W2$orderF)]
   
   b <- dose.5q31.W[, !is.na(sample.geno.W1$orderF)] 
   dose.5q31.W.F <- b[, order(sample.geno.W2$orderF)]
   
   b <- dose.15q26.W[, !is.na(sample.geno.W1$orderF)] 
   dose.15q26.W.F <- b[, order(sample.geno.W2$orderF)]
   
   b <- dose.1p13.W[, !is.na(sample.geno.W1$orderF)] 
   dose.1p13.W.F <- b[, order(sample.geno.W2$orderF)]
   
   b <- dose.15q24.W[, !is.na(sample.geno.W1$orderF)] 
   dose.15q24.W.F <- b[, order(sample.geno.W2$orderF)]
   
   b <- dose.1q41.W[, !is.na(sample.geno.W1$orderF)] 
   dose.1q41.W.F <- b[, order(sample.geno.W2$orderF)]
   
   ## merge back to gen expression file   
   pheno1.W1 <- merge(pheno1.W, sample.W.3, by.x="patientid", by.y="id_1", all=TRUE) 
   pheno1.W1 <- pheno1.W1[order(pheno1.W1$col2.x) , ] 
   pheno1.W2 <- pheno1.W1[!is.na(pheno1.W1$orderF) , ]
   
   tmp <- GeneExp.W[, !is.na(pheno1.W1$orderF)] 
   GeneExp.W.F <- tmp[, order(pheno1.W2$orderF)]
   
   ## merge back to CNV file      
   sample.CNV.W1 <- merge(sample.CNV.W, sample.W.3, by.x="ID.W", by.y="id_1", all=TRUE) 
   sample.CNV.W1 <- sample.CNV.W1[order(sample.CNV.W1$order3.x) , ] 
   sample.CNV.W2 <- sample.CNV.W1[!is.na(sample.CNV.W1$orderF) , ]
   
   CNV.W2 <- CNV.W.data4[, 4:ncol(CNV.W.data4)]
   
   tmp <- CNV.W2[, !is.na(sample.CNV.W1$orderF)] 
   CNV.W.F <- tmp[, order(sample.CNV.W2$orderF)]
   rownames(CNV.W.F) <- CNV.W.data4$GeneSymbol
   
   
### section [5] eQTL linear regression
## batch effect correction, with and without subtype adjustment
   sample.3 <- rbind(sample.B.3, sample.W.3)
   GeneExp.F <- cbind(GeneExp.B.F, GeneExp.W.F) 
   CNV.F <- as.matrix(cbind(CNV.B.F, CNV.W.F))
   
   dim(CNV.F)
   dim(GeneExp.F)
   dim(sample.3)
   
   sample.4 <- cbind(sample.3[, c("id_1", "age_at_diagnosis", "subtype_pam50", "batch")],
         data.frame(lapply(sample.3[, c("prop_yri", "prop_asn")], function(x) as.numeric(as.character(x)))))
   
   residual1 = NULL  
   residual2 = NULL     
   for (i in 1:nrow(GeneExp.F)){
      Data <- cbind.data.frame(exprs1=GeneExp.F[i, ], CNV1=CNV.F[i, ], sample.4)  
      
      lm.exp1 <- lm(exprs1 ~ CNV1 + age_at_diagnosis + as.factor(batch), data=Data) 
      res1 <- residuals(lm.exp1)
      residual1 <- rbind(residual1, res1)
      
      lm.exp2 <- lm(exprs1 ~ CNV1 + age_at_diagnosis + as.factor(batch) + as.factor(subtype_pam50), data=Data) 
      res2 <- residuals(lm.exp2)
      residual2 <- rbind(residual2, res2)
   }
   rownames(residual1) <- rownames(CNV.F)
   rownames(residual2) <- rownames(CNV.F)

      
## eQTL analysis
   sample.B.4 <- sample.4[1:140,] 
   sample.W.4 <- sample.4[141:883,]    
   residual1.B <- residual1[, 1:140]
   residual2.B <- residual2[, 1:140]
   residual1.W <- residual1[, 141:883]
   residual2.W <- residual2[, 141:883]
   
   loci <- unique(gene.anno$loci)
   dataName1 <- paste("result", "B", loci, sep=".")
   dataName2 <- paste("result", "W", loci, sep=".")
   index <- cbind.data.frame(loci, dataName1, dataName2, stringsAsFactors = FALSE)
   index$k <- 1:6
   index
   
   ## Blacks 
   for (k in 1:nrow(index)) { 
      
      express1 <- residual1.B[gene.anno$loci==index[k, "loci"], ]
      express2 <- residual2.B[gene.anno$loci==index[k, "loci"], ]
      
      exprs.anno <- gene.anno[gene.anno$loci==index[k, "loci"], ]

      if(index[k, "loci"]=="7q11") {GENO <- dose.7q11.B.F }
      if(index[k, "loci"]=="5q31") {GENO <- dose.5q31.B.F }
      if(index[k, "loci"]=="15q26") {GENO <- dose.15q26.B.F }
      if(index[k, "loci"]=="1p13") {GENO <- dose.1p13.B.F }
      if(index[k, "loci"]=="15q24") {GENO <- dose.15q24.B.F }
      if(index[k, "loci"]=="1q41") {GENO <- dose.1q41.B.F } 
      
      print(paste("eQTL ANALYSIS FOR", index[k, "loci"], "...", sep=" "))
      resultB1 <- NULL 
      resultB2 <- NULL 
      for (i in 1:nrow(express1)){
         exprs1 <- express1[i, ]
         exprs2 <- express2[i, ]
         exprs1.anno <- exprs.anno[i, ]

         resultA1 <- NULL 
         resultA2 <- NULL 
         for (j in 1:nrow(GENO)){
            geno1 <- GENO[j,]
            SNP1 <- rownames(GENO)[j]
            Data1 <- cbind(exprs1, geno1, sample.B.4)  
            Data2 <- cbind(exprs2, geno1, sample.B.4)  
            
            lm.1 <- lm(exprs1 ~ geno1 + prop_yri + prop_asn, data=Data1) 
            a <- cbind.data.frame(index[k, "loci"], exprs1.anno$gene_name, SNP1, t(summary(lm.1)$coefficients[2,]), model=1)
            resultA1 <- rbind(resultA1, a)
            
            lm.2 <- lm(exprs2 ~ geno1 + prop_yri + prop_asn, data=Data2) 
            b <- cbind.data.frame(index[k, "loci"], exprs1.anno$gene_name, SNP1, t(summary(lm.2)$coefficients[2,]), model=2)
            resultA2 <- rbind(resultA2, b)
            
         } ##j
         resultB1 <- rbind(resultB1, resultA1) 
         resultB2 <- rbind(resultB2, resultA2) 
         print(paste("Finish", nrow(GENO), "regressions for", exprs1.anno$gene_name, sep=" "))
      } ##i
      
      resultB <- rbind(resultB1, resultB2)
      colnames(resultB) <- c("loci", "geneSymbol", "SNP", "beta", "SE", "T_value", "P_value", "model") 
      resultC <- resultB[resultB$P_value<0.05, c(1:5,7:8)]

      outfile1 <- paste("/gpfs/data/huo-lab/ROOT/Tunde/eQTL/new_result_eQTL_black_", index[k, "loci"], ".txt.gz", sep="")
      outfile2 <- paste("/gpfs/data/huo-lab/ROOT/Tunde/eQTL/new_result005_eQTL_black_", index[k, "loci"], ".txt", sep="")      
      write.table(resultB, gzfile(outfile1), quote =FALSE, row.names=FALSE, col.names=TRUE) 
      write.table(resultC, outfile2, quote =FALSE, row.names=FALSE, col.names=TRUE)       

   } ##k
   
## Whites 
   for (k in 1:nrow(index)) { 
      
      express1 <- residual1.W[gene.anno$loci==index[k, "loci"], ]
      express2 <- residual2.W[gene.anno$loci==index[k, "loci"], ]
      
      exprs.anno <- gene.anno[gene.anno$loci==index[k, "loci"], ]
      
      if(index[k, "loci"]=="7q11") {GENO <- dose.7q11.W.F }
      if(index[k, "loci"]=="5q31") {GENO <- dose.5q31.W.F }
      if(index[k, "loci"]=="15q26") {GENO <- dose.15q26.W.F }
      if(index[k, "loci"]=="1p13") {GENO <- dose.1p13.W.F }
      if(index[k, "loci"]=="15q24") {GENO <- dose.15q24.W.F }
      if(index[k, "loci"]=="1q41") {GENO <- dose.1q41.W.F } 
      
      print(paste("eQTL ANALYSIS FOR", index[k, "loci"], "...", sep=" "))
      resultB1 <- NULL 
      resultB2 <- NULL 
      for (i in 1:nrow(express1)){
         exprs1 <- express1[i, ]
         exprs2 <- express2[i, ]
         exprs1.anno <- exprs.anno[i, ]
         
         resultA1 <- NULL 
         resultA2 <- NULL 
         for (j in 1:nrow(GENO)){
            geno1 <- GENO[j,]
            SNP1 <- rownames(GENO)[j]
            Data1 <- cbind(exprs1, geno1, sample.W.4)  
            Data2 <- cbind(exprs2, geno1, sample.W.4)  
            
            lm.1 <- lm(exprs1 ~ geno1 + prop_yri + prop_asn, data=Data1) 
            a <- cbind.data.frame(index[k, "loci"], exprs1.anno$gene_name, SNP1, t(summary(lm.1)$coefficients[2,]), model=1)
            resultA1 <- rbind(resultA1, a)
            
            lm.2 <- lm(exprs2 ~ geno1 + prop_yri + prop_asn, data=Data2) 
            b <- cbind.data.frame(index[k, "loci"], exprs1.anno$gene_name, SNP1, t(summary(lm.2)$coefficients[2,]), model=2)
            resultA2 <- rbind(resultA2, b)
            
         } ##j
         resultB1 <- rbind(resultB1, resultA1) 
         resultB2 <- rbind(resultB2, resultA2) 
         print(paste("Finish", nrow(GENO), "regressions for", exprs1.anno$gene_name, sep=" "))
      } ##i
      
      resultB <- rbind(resultB1, resultB2)
      colnames(resultB) <- c("loci", "geneSymbol", "SNP", "beta", "SE", "T_value", "P_value", "model") 
      resultC <- resultB[resultB$P_value<0.05, c(1:5,7:8)]
      
      outfile1 <- paste("/gpfs/data/huo-lab/ROOT/Tunde/eQTL/new_result_eQTL_white_", index[k, "loci"], ".txt.gz", sep="")
      outfile2 <- paste("/gpfs/data/huo-lab/ROOT/Tunde/eQTL/new_result005_eQTL_white_", index[k, "loci"], ".txt", sep="")      
      write.table(resultB, gzfile(outfile1), quote =FALSE, row.names=FALSE, col.names=TRUE) 
      write.table(resultC, outfile2, quote =FALSE, row.names=FALSE, col.names=TRUE)       

   } ##k

   
   gene.anno.out <- gene.anno[, -9]
   geno.1p13.left$V1 <- "1p13"
   geno.1q41.left$V1 <- "1q41"
   geno.5q31.left$V1 <- "5q31"
   geno.7q11.left$V1 <- "7q11"
   geno.15q24.left$V1 <- "15q24"
   geno.15q26.left$V1 <- "15q26"
   geno.left <- rbind(geno.1p13.left, geno.1q41.left, geno.5q31.left, geno.7q11.left, geno.15q24.left, geno.15q26.left)
   
   outfile3 <- "/gpfs/data/huo-lab/ROOT/Tunde/eQTL/gene.226.anno.txt"
   outfile4 <- "/gpfs/data/huo-lab/ROOT/Tunde/eQTL/SNP.6loci.anno.txt"
   write.table(gene.anno.out, outfile3, quote =FALSE, row.names=FALSE, col.names=TRUE) 
   write.table(geno.left, outfile4, quote =FALSE, row.names=FALSE, col.names=TRUE)    
                          
      
sink() 
