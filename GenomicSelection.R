#install.packages("rrBLUP")
#install.packages("multtest")
#install.packages("gplots")
#install.packages("LDheatmap")
#install.packages("genetics")
#install.packages("EMMREML")
#install.packages("compiler")

library(rrBLUP)
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) 
library("scatterplot3d")

source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

nam3_hmp <- read.delim("C:/Users/yismail/Desktop/C.C.A/All.Chromosomes.YI.2019.hmp.txt", head=F) 

myG <- nam3_hmp
myGAPIT <- GAPIT(G=myG, output.numerical=TRUE, SNP.impute = "Middle", SNP.MAF=0.03)
dat=myGAPIT$GD

nan=dat[,1]
nan <- as.vector(as.matrix(nan))
dat <- dat[,2:ncol(dat)]
dat=as.matrix(dat)
#dat[which(dat==2)]=1
rownames(dat)=nan
##
nam <- data.frame(nan)
names(nam) <- "Taxa"

df2 <- as.matrix(dat)
df2 <- df2-1
df2[1:6,1:6]

pheno <- read.csv("C:/Users/yismail/Desktop/C.C.A/Nirdata.csv", header=T)
head(pheno)

# Create Taxa column and remove flowcell information from accession names
pheno$Taxa <- sub("\\:.*", "", pheno$Taxa)
head(pheno)
nrow(pheno)

phenames <- colnames(pheno[,-1])

nam <- data.frame(rownames(df2))
names(nam) <- "Taxa"
Pheno <- merge(pheno, nam, by="Taxa")
df2.copy <- df2 
df2 <- df2[match(Pheno$Taxa, rownames(df2)),]

# There is no need for imputation since our GBS data had been previously imputed
#impute=A.mat(df2,max.missing=0.5,impute.method="mean",return.imputed=T)

cycles = 50
k = 5
trt_accuracy.df.null=NULL

for(p in 1:length(phenames)){
  trt_accuracy=matrix(nrow=cycles, ncol=k, NA)
  
  CrossPheno <- data.frame(Pheno[,1], Pheno[,phenames[p]])
  colnames(CrossPheno)[2] <- phenames[p]
  colnames(CrossPheno)[1] <- "Taxa"
  
  
  for(j in 1:cycles){
    CrossPheno$id <- sample(1:k, nrow(CrossPheno), replace = TRUE)
    list <- 1:k
    
    for (i in 1:k){
      
      # remove rows with id i from dataframe to create training set
      # select rows with id i to create test set
      pheno_trainingset <- subset(CrossPheno, id %in% list[-i])
      Geno_trainingset <- df2[match(pheno_trainingset$Taxa, rownames(df2)),]
      
      pheno_testset <- subset(CrossPheno, id %in% c(i))
      Geno_testset <- df2[match(pheno_testset$Taxa, rownames(df2)),]
      
      # run RRBLUP
      trt=(pheno_trainingset[,2])
      trt_answer<-mixed.solve(trt, Z=Geno_trainingset, K=NULL, SE = FALSE, return.Hinv=FALSE)
      trt_u = trt_answer$u
      e_trt = as.matrix(trt_u)
      pred_trt_valid =  Geno_testset %*% e_trt
      pred_trt=(pred_trt_valid[,1])+ trt_answer$beta
      
      trt_testset = pheno_testset[,2]
      trt_accuracy[j,i] <- cor(pred_trt, trt_testset, use="complete" )
      
      cat(".")
      
    }
    cat(">")
  }
  trt_accuracy.df <- data.frame(trt_accuracy)
  trt_accuracy.df$Trait <- rep(phenames[p], nrow(trt_accuracy.df))
  trt_accuracy.df.null <- rbind(trt_accuracy.df.null, trt_accuracy.df)
  trt_accuracy <- NULL
  trt_accuracy.df <- NULL
}

write.table(trt_accuracy.df.null, "Maize.AMES.GrainQuality.11112019.csv", sep=",", quote=F, row.names=F, col.names=T)

