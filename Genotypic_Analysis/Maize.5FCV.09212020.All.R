rm(list=ls())
setwd("/homes/omo/public_html/vGWAS/Maize.V2/GS_New/")
date="10262020" 
#setwd("")
library(rrBLUP)
library(kernlab)
library(e1071)
library(randomForest)
library(BGLR)
library(sommer)
source("http://people.beocat.ksu.edu/~omo/vGWAS/Maize.V2/GS_New/Maize.GS.RRBLUP.RKHS.SVR.RF.Algorithm.09212020.R")
source("https://raw.githubusercontent.com/ekfchan/evachan.org-Rscripts/master/rscripts/geno_to_allelecnt.R")
source("http://people.beocat.ksu.edu/~omo/Collaborations/Sorghum/CoincidenceIndex.R")

GD <- read.delim("http://people.beocat.ksu.edu/~omo/vGWAS/Maize.V2/GS_New/Maize.Numerical.Genotypes.953Inds.264428SNPs.txt", header=T, stringsAsFactors = F, skip = 1)
#GD <- read.delim("/Volumes/Seagate Expansion Drive/Mac Files/PostDoc/Yasser/Maize.Panzea/GS_New/Maize.Numerical.Genotypes.953Inds.264428SNPs.txt", header=T, stringsAsFactors = F, skip = 1)

GD[1:6,1:6]
#dat=myGAPIT$GD
dat=GD
nan=dat[,1]
nan <- as.vector(as.matrix(nan))
dat <- dat[,2:ncol(dat)]
dat=as.matrix(dat)
#dat=dat-1
rownames(dat)=nan
dat[1:6,1:6]
dat[dat==1] <- 2
dat[dat==0.5] <- 1
dat[1:6,1:6]

df.nam.geno <- data.frame(rownames(dat))
names(df.nam.geno) <- "Taxa"

df2 <- as.matrix(dat)
#df2-1
df2[1:6,1:6]
any(is.na(df2))
sum(is.na(df2))

dim(df.nam.geno)
df.nam.geno2 <- as.data.frame(df.nam.geno[!duplicated(df.nam.geno$Taxa),])
dim(df.nam.geno2)
head(df.nam.geno2)
colnames(df.nam.geno2) <- "Taxa"
head(df.nam.geno2)

phen.taxa <- read.csv("http://people.beocat.ksu.edu/~omo/vGWAS/Nirdata.csv", header=T)
dim(phen.taxa)
phen.taxa2 <- phen.taxa[!duplicated(phen.taxa$Taxa),]
dim(phen.taxa2)
head(phen.taxa2)

# PI559935, Ames26133

colnames(phen.taxa2)[1] <- "Taxa"
phen.taxa2$Taxa <- as.character(phen.taxa2$Taxa)
head(phen.taxa2)

com.taxa <- merge(phen.taxa2, df.nam.geno2)
dim(com.taxa)
head(com.taxa)
com.taxa$Taxa <- as.character(com.taxa$Taxa)

com.taxa.df <- data.frame(com.taxa$Taxa)
colnames(com.taxa.df) <- "Taxa"

head(com.taxa)
names(com.taxa)
com.taxa2 <- com.taxa[,-3]
head(com.taxa2)
str(com.taxa2)
phenames <- as.character(colnames(com.taxa2[,-1]))
phenames <- phenames

df3 <- df2[match(com.taxa2$Taxa, rownames(df2)),]
any(is.na(df3))
k = 5
cycles=10

number.of.folds=k

species <- "Maize"




diploid.polyploid.5cfv(phenames=phenames, m_train.pheno=com.taxa2, diploid.Geno=df3, species=species, k=k, cycles = cycles,
                                  number.of.folds=number.of.folds, date=date)



phenames=phenames; m_train.pheno=com.taxa2; diploid.Geno=df3; species=species; k=k; cycles = cycles; number.of.folds=number.of.folds; date=date
any(is.na(df3))

