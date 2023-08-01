rm(list=ls())
#setwd("~/Google Drive/Post Doc/Collaborative Publications/Cowpea/Analysis/Results")
setwd("/homes/omo/public_html/vGWAS/Brians_Data/MLMM.2p/")

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

myG <- read.delim("http://people.beocat.ksu.edu/~omo/vGWAS/Brians_Data/Yasser.Maize.LDImpute.v1.866Inds.303489SNPs.2pMAF.hmp.txt", header=T)
myGAPIT <- GAPIT(G=myG, output.numerical=TRUE)
dat= myGAPIT$GD #read.delim("/homes/omo/public_html/vGWAS/", header = T, stringsAsFactors = F)#
dat <- read.delim("/homes/omo/public_html/vGWAS/Maize.V2/GAPIT.Genotype.Numerical.txt", header=T, stringsAsFactors=F)
dat[1:6,1:6]
nan=dat[,1]
nan <- as.vector(as.matrix(nan))
dat <- dat[,2:ncol(dat)]
dat=as.matrix(dat)

rownames(dat)=nan
##
nam <- data.frame(nan)
names(nam) <- "Taxa"

dat[1:6,1:6]

##
nam <- data.frame(nan)
names(nam) <- "Taxa"
nam$Taxa <- as.character(nam$Taxa)
nam <- nam[!duplicated(nam$Taxa),]
nam <- data.frame(nam)
names(nam) <- "Taxa"

df2 <- as.matrix(dat)
#df2 <- df2-1
df2[1:6,1:6]

pheno <- read.csv("http://people.beocat.ksu.edu/~omo/vGWAS/Nirdata.csv", header=T)
head(pheno)
pheno <- pheno[!duplicated(pheno$Taxa),]
# Create Taxa column and remove flowcell information from accession names


Phenotypes <- merge(pheno, nam)


phenames <- colnames(Phenotypes[,-c(1,3)]) # Density and Fibre have no significant associations


df2 <- df2[match(Phenotypes$Taxa, rownames(df2)),]
  dim(df2)

for(p in 1:length(phenames)){
      #df2 <- df2[!duplicated(rownames(df2)),]
      JL_RES <- read.csv(paste("/homes/omo/public_html/vGWAS/Maize.V2/MLMM/GWAS_mlmm_results.Konly",phenames[p],".csv", sep=""), header=T, stringsAsFactors = F)
      head(JL_RES)
      JL_RES$Trait <- rep(phenames[p], nrow(JL_RES))
      JL_RES$SNP <- as.character(JL_RES$SNP)
      Bonf <- 0.05/nrow(JL_RES)
      JL_RES <- JL_RES[which(JL_RES$pval < Bonf),]
      JL_RES
      nrW <- nrow(JL_RES)
      if(nrW==0){print(paste("--- No Significant QTL Identified For Trait: ", phenames[p], " ---", sep = "")) ;next}else{
      #phenames <- names(com.tax[,-1])
      #for (l in 1:length(phenames)){
        print(paste("-------------- Trait being analysed: ", phenames[p], "!!!!!!!!!!!---------------", sep = ""))
        
        ExplVar200Best <- JL_RES#[which(JL_RES$Trait==phenames[p]),]
        snp.vec <- ExplVar200Best$SNP
        #bSNP<-df2[,as.character(ExplVar200Best$SNP)]
        
        if(length(snp.vec)==1){# If there is only one vQTL
          bSNP<-df2[,as.character(ExplVar200Best$SNP)]
          bSNP <- data.frame(bSNP)
          colnames(bSNP) <- snp.vec
          bSNP[,1] <- as.character(bSNP[,1])
          
        }else{bSNP<-df2[,as.character(ExplVar200Best$SNP)]}
        
        
        phdata <- data.frame(Phenotypes[,1], Phenotypes[,phenames[p]])
        colnames(phdata)[2] <- phenames[p]
        colnames(phdata)[1] <- "Taxa"
        #sP<-as.data.frame(phdata[,phenames[l]])
        sP<-phdata
        sP <- sP[!duplicated(sP$Taxa),]
        rownames(sP) <- sP$Taxa
        
        #bSNP2 <- bSNP[!duplicated(names(bSNP))]
        #bSNP2 <- bSNP[!duplicated(bSNP)]
        da<-as.data.frame(cbind(sP, bSNP))
        
        trait_QTL_Pheno <- da
        write.table(t(data.frame(c("Trait", "QTL", "Additive Effect", "PVE"))), paste("Maize.QTL.Effects_", phenames[p],"_QTL",".txt", sep=""), sep="\t", append=T, quote=F, row.names=F, col.names=F)
        #APV is the Among population variance in accordance to Wurschum et al. 2011 Heredity
        for(i in 3:ncol(trait_QTL_Pheno)){
          
          snp <- colnames(trait_QTL_Pheno)[i]  
          print(paste("-------------- Trait being analysed: ", phenames[p], "SNP: ", snp, "!!!!!!!!!!!---------------", sep = ""))
          
          trait_QTL_Pheno_2 <- trait_QTL_Pheno[,c(1,2,i)]
          
          AA_class <- trait_QTL_Pheno[which(trait_QTL_Pheno[,i]==2),]
          AA <- mean(AA_class[,2], na.rm=T)
          BB_class <- trait_QTL_Pheno[which(trait_QTL_Pheno[,i]==0),]
          BB <- mean(BB_class[,2], na.rm=T)
          QTL_effect <- (AA-BB)/2
          
          #formula.single <- as.formula(paste("Cd_comb ~ ",paste(as.character(topSNP$SNP), collapse=" + "), sep=" "))
          trait_QTL_Pheno_2$QTL <- trait_QTL_Pheno_2[,3]
          #QTL <- colnames(trait_QTL_Pheno_2[3])
          fin.anova <- lm(trait_QTL_Pheno_2[,phenames[p]] ~  QTL, data=trait_QTL_Pheno_2, na.action = na.omit)
          fin.sum <- summary(fin.anova)
          
          QVar <- round((fin.sum$adj.r.squared)*100, digits=2)#Phenotypes.FT[,phenames[l]]
          
          print(paste("-------------- PVE For SNP: ", snp, "; Trait: ", phenames[p], " == ", QVar, "%  !!!!!!!!!!!---------------", sep = ""))
          write.table(t(data.frame(c(phenames[p], colnames(trait_QTL_Pheno[i]), round(abs(QTL_effect[1]), 3), QVar[1]))), paste("Maize.QTL.Effects_", phenames[p],"_QTL",".txt", sep=""), sep="\t", append=T, quote=F, row.names=F, col.names=F)
        }
        
        
      }

}


Ash.PVE <- read.delim("/homes/omo/public_html/vGWAS/Maize.V2/MLMM/Maize.QTL.Effects_ASH_QTL.txt", header=T, stringsAsFactors = F)
Den.PVE <- read.delim("/homes/omo/public_html/vGWAS/Maize.V2/MLMM/Maize.QTL.Effects_DEN_QTL.txt", header=T, stringsAsFactors = F)
Fib.PVE <- read.delim("/homes/omo/public_html/vGWAS/Maize.V2/MLMM/Maize.QTL.Effects_FIB_QTL.txt", header=T, stringsAsFactors = F)
Pro.PVE <- read.delim("/homes/omo/public_html/vGWAS/Maize.V2/MLMM/Maize.QTL.Effects_PRO_QTL.txt", header=T, stringsAsFactors = F)
Sta.PVE <- read.delim("/homes/omo/public_html/vGWAS/Maize.V2/MLMM/Maize.QTL.Effects_STA_QTL.txt", header=T, stringsAsFactors = F)
QTLs <- rbind(Ash.PVE, Den.PVE, Fib.PVE, Pro.PVE, Sta.PVE) #Oil.PVE,
QTLs$QTL <- as.character(QTLs$QTL)
df3 <- df2[,match(QTLs$QTL, colnames(df2))]
df3[1:6,1:5]
genot <- t(df3)
genot[1:5,1:5]
rownames(genot)
genot2 <- genot[!duplicated(rownames(genot)),]
genot2[1:5,1:5]
source("https://raw.githubusercontent.com/ekfchan/evachan.org-Rscripts/master/rscripts/calc_snp_stats.R")
#geno = p x n
#class(geno)

snp.summary <- calc_snp_stats(geno=genot2)
snp.summary$QTL <- rownames(snp.summary)
snp.summary2 <- snp.summary[,c(14,6)]
QTLs2 <- merge(QTLs, snp.summary2)
write.table(QTLs2, "/homes/omo/public_html/vGWAS/Maize.V2/MLMM/QTL.Summary.Maize.08282020.csv", sep=",", quote=F, row.names = F, col.names = T)



#####
#### vQTL Analysis
setwd("/homes/omo/public_html/vGWAS/Maize.V2/vGWAS/")
head(Phenotypes)
phenames <- colnames(Phenotypes[,-c(1,3)]) # Density and Fibre have no significant associations


df2 <- df2[match(Phenotypes$Taxa, rownames(df2)),]
dim(df2)

for(p in 1:length(phenames)){
  #df2 <- df2[!duplicated(rownames(df2)),]
  JL_RES <- read.csv(paste("/homes/omo/public_html/vGWAS/Maize.V2/vGWAS/vQTLs.Effect",phenames[p],"csv", sep="."), header=T, stringsAsFactors = F)
  head(JL_RES)
  colnames(JL_RES)[3] <- "SNP"
  JL_RES$Trait <- rep(phenames[p], nrow(JL_RES))
  JL_RES$SNP <- as.character(JL_RES$SNP)

  print(paste("-------------- Trait being analysed: ", phenames[p], "!!!!!!!!!!!---------------", sep = ""))
  
  ExplVar200Best <- JL_RES#[which(JL_RES$Trait==phenames[p]),]
  
  bSNP<-df2[,as.character(ExplVar200Best$SNP)]
  
  phdata <- data.frame(Phenotypes[,1], Phenotypes[,phenames[p]])
  colnames(phdata)[2] <- phenames[p]
  colnames(phdata)[1] <- "Taxa"
  #sP<-as.data.frame(phdata[,phenames[l]])
  sP<-phdata
  sP <- sP[!duplicated(sP$Taxa),]
  rownames(sP) <- sP$Taxa
  
  #bSNP2 <- bSNP[!duplicated(names(bSNP))]
  #bSNP2 <- bSNP[!duplicated(bSNP)]
  da<-as.data.frame(cbind(sP, bSNP))
  
  trait_QTL_Pheno <- da
  head(trait_QTL_Pheno)
  write.table(t(data.frame(c("Trait", "QTL", "Additive Effect", "PVE", "Vm", "Vv"))), paste("Maize.vQTL.Effects_", phenames[p],"_QTL",".txt", sep=""), sep="\t", append=T, quote=F, row.names=F, col.names=F)
  #APV is the Among population variance in accordance to Wurschum et al. 2011 Heredity
  for(i in 3:ncol(trait_QTL_Pheno)){
    
    snp <- colnames(trait_QTL_Pheno)[i]  
    print(paste("-------------- Trait being analysed: ", phenames[p], "SNP: ", snp, "!!!!!!!!!!!---------------", sep = ""))
    
    trait_QTL_Pheno_2 <- trait_QTL_Pheno[,c(1,2,i)]
    
    AA_class <- trait_QTL_Pheno[which(trait_QTL_Pheno[,i]==2),]
    AA <- mean(AA_class[,2], na.rm=T)
    BB_class <- trait_QTL_Pheno[which(trait_QTL_Pheno[,i]==0),]
    BB <- mean(BB_class[,2], na.rm=T)
    QTL_effect <- (AA-BB)/2
    
    #formula.single <- as.formula(paste("Cd_comb ~ ",paste(as.character(topSNP$SNP), collapse=" + "), sep=" "))
    trait_QTL_Pheno_2$QTL <- trait_QTL_Pheno_2[,3]
    #QTL <- colnames(trait_QTL_Pheno_2[3])
    fin.anova <- lm(trait_QTL_Pheno_2[,phenames[p]] ~  QTL, data=trait_QTL_Pheno_2, na.action = na.omit)
    fin.sum <- summary(fin.anova)
    
    QVar <- round((fin.sum$adj.r.squared)*100, digits=2)#Phenotypes.FT[,phenames[l]]
    
    print(paste("-------------- PVE For SNP: ", snp, "; Trait: ", phenames[p], " == ", QVar, "%  !!!!!!!!!!!---------------", sep = ""))
    write.table(t(data.frame(c(phenames[p], colnames(trait_QTL_Pheno[i]), round(abs(QTL_effect[1]), 3), QVar[1], round(abs(JL_RES[(i-2),1]),5), round(abs(JL_RES[(i-2),2]),5)))), paste("Maize.vQTL.Effects_", phenames[p],"_QTL",".txt", sep=""), sep="\t", append=T, quote=F, row.names=F, col.names=F)

  }
  
  
  
}

rm(df3)
Ash.PVE <- read.delim("/homes/omo/public_html/vGWAS/Maize.V2/vGWAS/Maize.vQTL.Effects_ASH_QTL.txt", header=T, stringsAsFactors = F)
Ash.PVE$Trait <- rep("ASH", nrow(Ash.PVE))
Den.PVE <- read.delim("/homes/omo/public_html/vGWAS/Maize.V2/vGWAS/Maize.vQTL.Effects_DEN_QTL.txt", header=T, stringsAsFactors = F)
Den.PVE$Trait <- rep("DEN", nrow(Den.PVE))
Pro.PVE <- read.delim("/homes/omo/public_html/vGWAS/Maize.V2/vGWAS/Maize.vQTL.Effects_PRO_QTL.txt", header=T, stringsAsFactors = F)
Pro.PVE$Trait <- rep("PRO", nrow(Pro.PVE))
Sta.PVE <- read.delim("/homes/omo/public_html/vGWAS/Maize.V2/vGWAS/Maize.vQTL.Effects_STA_QTL.txt", header=T, stringsAsFactors = F)
Sta.PVE$Trait <- rep("STA", nrow(Sta.PVE))
Fib.PVE <- read.delim("/homes/omo/public_html/vGWAS/Maize.V2/vGWAS/Maize.vQTL.Effects_FIB_QTL.txt", header=T, stringsAsFactors = F)
Fib.PVE$Trait <- rep("FIB", nrow(Fib.PVE))
QTLs <- rbind(Ash.PVE, Den.PVE, Pro.PVE, Sta.PVE, Fib.PVE)#, Oil.PVE
head(QTLs)
colnames(QTLs)[2] <- "QTL"
QTLs$QTL <- as.character(QTLs$QTL)
df3 <- df2[,match(QTLs$QTL, colnames(df2))]
df3[1:6,1:6]
genot <- t(df3)
genot[1:6,1:6]
rownames(genot)
genot2 <- genot[!duplicated(rownames(genot)),]
genot2[1:6,1:6]
source("https://raw.githubusercontent.com/ekfchan/evachan.org-Rscripts/master/rscripts/calc_snp_stats.R")
#geno = p x n
#class(geno)

snp.summary <- calc_snp_stats(geno=genot2)
snp.summary$QTL <- rownames(snp.summary)
snp.summary2 <- snp.summary[,c(14,6)]
QTLs2 <- merge(QTLs, snp.summary2)
head(QTLs2)
write.table(QTLs2, "/homes/omo/public_html/vGWAS/Maize.V2/vGWAS/vQTL.Summary.Maize.08282020.csv", sep=",", quote=F, row.names = F, col.names = T)



GWAS.loci <- read.csv("/Users/omo/Desktop/PostDoc/Yasser/NewMLMM/QTL.Summary.Maize.06072020.csv", header = T, stringsAsFactors = F)
head(GWAS.loci)
GWAS.loci$Model <- rep("GWAS", nrow(GWAS.loci))
vGWAS.loci <- read.csv("/Users/omo/Desktop/PostDoc/Yasser/NewMLMM/vQTL.Summary.Maize.06072020.csv", header = T, stringsAsFactors = F)
vGWAS.loci$Model <- rep("vGWAS", nrow(vGWAS.loci))
head(vGWAS.loci)
Loci.Associatns <-rbind(GWAS.loci, vGWAS.loci)








