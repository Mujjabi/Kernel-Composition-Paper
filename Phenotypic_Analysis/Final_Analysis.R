#This is the final workflow for the phenotypic analysis of kernel compositional traits
# I conducted an analysis of variance between heterotic groups since this wasnt possible between inbreds due to lack of reps
#I used a simple linear model which assumes fixed factors, since we dont have random factors. Heterotic groups are  fixed groups thus considered fixed

Compiled <- read.csv("GrainQualityData_AmesPanel.csv", header = TRUE)

Compiled$Accession <- as.factor(Compiled$Accession)
Compiled$Origin <- as.factor(Compiled$Origin)
Compiled$Group <- as.factor(Compiled$Group)
Compiled$Starch <- as.numeric(Compiled$Starch)
Compiled$Oil <- as.numeric(Compiled$Oil)
Compiled$Protein <- as.numeric(Compiled$Protein)
Compiled$Fiber <- as.numeric(Compiled$Fiber)
Compiled$Density <- as.numeric(Compiled$Density)
Compiled$Ash <- as.numeric(Compiled$Ash)



## Testing for outliers
boxplot(Compiled[4:9], plot=TRUE)$out


# create detect outlier function

outliers <- function(x) {
  
  Q1 <- quantile(x, probs=.25) # calculate first quantile
  Q3 <- quantile(x, probs=.75) # calculate third quantile
  iqr = Q3-Q1           # calculate inter quartile range
  
  upper_limit = Q3 + (iqr*1.5)  #Calculate upper limit
  lower_limit = Q1 - (iqr*1.5)  #calculate lower limit
  
  x > upper_limit | x < lower_limit   # return true or false
}

remove_outliers <- function(Compiled, cols = names()) {  # for loop to traverse in columns vector
  for (col in cols) {
    Compiled <- Compiled[!outliers(Compiled[[col]]),] # remove observation if it satisfies outlier function
  }
  Compiled
}

#Remove outliers 
Compiled2 <- remove_outliers(Compiled, c(4:9))

outliers <- boxplot(Compiled[4:9], plot=TRUE)$out  # save the outliers in a vector
x<-Compiled
Compiled <- x[-which(x$Starch %in% outliers),] #Removing outlier for just one variable
Compiled <- Compiled[-which(Compiled$Protein %in% outliers),]
Compiled <- Compiled[-which(Compiled$Oil %in% outliers),]
Compiled <- Compiled[-which(Compiled$Fib %in% outliers),]
Compiled <- Compiled[-which(Compiled$Ash %in% outliers),]
Compiled <- Compiled[-which(Compiled$Density %in% outliers),]

boxplot(Compiled2[4:9], col=c("red", "yellow", "green","#999999", "#E69F00", "#56B4E9"), plot=TRUE)$out 

### Only model with Group can be used. We have reps within groups. 

modelprot <- aov(Protein ~ Group, Compiled2)
modelstr <- aov(Starch ~ Group, Compiled2)
modeloil <- aov(Oil ~ Group, Compiled2)
modelash <- aov(Ash ~ Group, Compiled2)
modeldensity <- aov(Density ~ Group, Compiled2)
modelfiber <- aov(Fib ~ Group, Compiled2)

anova(modelprot)
anova(modelstr)
anova(modeloil)
anova(modelfiber)
anova(modelash)
anova(modeldensity)

## Pairwise Mean Comparisons for each trait. 
## This is Tukeys HSD Test

library(agricolae)
out<-HSD.test(modelprot,"Group", group=TRUE)
print(out)

# Old version HSD.test()
df<-df.residual(modelprot)
MSerror<-deviance(modelprot)/df
with(Compiled2,HSD.test(Protein,Group,df,MSerror, group=TRUE,console=TRUE,
                          main="Protein Composition for Different Maize Groups "))


### spearman Rank Correlations between our values and grain quality values obtained by Hirsch et al, 2021. 
## Adjusted Blup values were extracted from Table 2 of the paper, and correlated to our estimated blups for each group 

#First, we estimate trait Blups for each Group using a mixed model
### Mixed Models 

library(lme4)
library(jtools)
library(lmerTest)
library(car)
install.packages("xlsx")
library(xlsx)
library(readr)

mod1 <- lmer(Protein ~ (1|Group), REML = TRUE, data= Compiled)
mod2 <- lmer(Starch ~ (1|Group), REML = TRUE, data= Compiled)
mod3 <- lmer(Oil ~ (1|Group), REML = TRUE, data= Compiled)
mod4 <- lmer(Ash ~ (1|Group), REML = TRUE, data= Compiled)
mod5 <- lmer(Fib ~ (1|Group), REML = TRUE, data= Compiled)
mod6 <- lmer(Density ~ (1|Group), REML = TRUE, data= Compiled)


## extract blups from each model
varComp<-as.data.frame(VarCorr(mod1,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
blupProt= coef(mod1)$Group

varComp<-as.data.frame(VarCorr(mod3,comp="vcov"))
blupOil = coef(mod3)$Group

varComp<-as.data.frame(VarCorr(mod2,comp="vcov"))
blupStr = coef(mod2)$Group

varComp<-as.data.frame(VarCorr(mod5,comp="vcov"))
blupfiber = coef(mod5)$Group

varComp<-as.data.frame(VarCorr(mod4,comp="vcov"))
blupash = coef(mod4)$Group

varComp<-as.data.frame(VarCorr(mod6,comp="vcov"))
blupdensity = coef(mod6)$Group

Ames_Blups <- cbind(blupProt,blupStr,blupOil,blupfiber,blupdensity,blupash)

write.csv(Ames_Blups, "Ames_Blups.csv")


# Test for normality Assumptioms
# qq plot
par(mfrow=c(3,2))
qqPlot(residuals(mod1), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Starch Content") 
qqPlot(residuals(mod2), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Protein Content", id = FALSE) 
qqPlot(residuals(mod3), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Oil Content", id = FALSE) 
qqPlot(residuals(mod4), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Ash Content", id = FALSE) 
qqPlot(residuals(mod5), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Fiber Content", id = FALSE) 
qqPlot(residuals(mod6), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Kernel Density", id = FALSE) 

#These Blups were combined with the Hirsch et al, 2021 Blups from Table 2 to for "TableValues" dataset used in the next step
#The dataset has the Group, Trait, Ames_Blups for each trait and Wisc_Blups for each trait. 
#Group performance/rankig for a particular trait across the studies was evaluated using a spearman ranking test below. 

Data1  <- read.csv("TableValue.csv", header = TRUE)

library(dplyr) 
Corr1 <- Data1%>% 
  group_by (Trait) %>% 
  summarise(cor=cor(Ames_Blups,Wisc_Blups, method = "spearman"))

write.csv(Corr1, "spearmanCorr.csv")

### Pearson correlation using mean values of the 275 genotypes that where phenotyped in both studies

Wis <- read.csv("Wis-Rawdata.csv", header = TRUE)  #Subset of Hirsch et al, 2021 raw-dataset with only 275 common, from all tested environment 

Wis$Accession <- as.factor(Wis$Accession)
Wis$Genotype <- as.factor(Wis$Genotype)
Wis$Group <- as.factor(Wis$Group)
Wis$Starch <- as.numeric(Wis$Starch)
Wis$Oil <- as.numeric(Wis$Oil)
Wis$Prot <- as.numeric(Wis$Prot)
Wis$Fiber <- as.numeric(Wis$Fiber)
Wis$Ash <- as.numeric(Wis$Ash)

library(dplyr)
WisMeans <- Wis %>%
  group_by(Accession) %>%
  summarise(across(Prot_W:Star_W, mean)) #Means of the 275 genotypes 
write.csv (WisMeans, "WisMeans.csv") 

##Combining the 2 datsets from the 2 studies 
#Combined Dataset with grain quality trait values from both studies

WisAmesCombined <- merge(Compiled, WisMeans, by = 'Accession') #merge function joins the 2 datasets by the common accession 
write.csv (WisAmesCombined, "Combined2.csv") 

## Pearson Correlation matrix of the observed and literature values for the 275 common genotypes. 

library(Hmisc) #Create the function to generate the correlation matrix
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )}
res3<-rcorr(as.matrix(WisAmesCombined[,4:14]),type = "pearson")
round(res3$r, 2)
round(res3$P, 2)
Corr3 <- flattenCorrMatrix(res3$r, res3$P)
write.csv (Corr2, "PearsonCorr.csv")


#Plotting Correlations

library("ggpubr")
Starch <- ggscatter(WisAmesCombined, x = "Starch", y = "Star_W", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          title = "Starch",
          xlab = "Observed [%]", ylab = "Literature [%]")

Protein <- ggscatter(WisAmesCombined, x = "Protein", y = "Prot_W", 
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson",
                    title = "Protein",
                    xlab = "Observed [%]", ylab = "Literature  [%]")

Oil <- ggscatter(WisAmesCombined, x = "Oil", y = "Oil_W", 
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson",
                    title = "Oil",
                 xlab = "Observed [%]", ylab = "Literature  [%]")

Fiber <- ggscatter(WisAmesCombined, x = "Fib", y = "Fib_W", 
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson",
                    title = "Fiber",
                    xlab = "Observed [%]", ylab = "Literature  [%]")
                  

Ash <- ggscatter(WisAmesCombined, x = "Ash", y = "Ash_W", 
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "pearson",
                    title = "Ash",
                    xlab = "Observed [%]", ylab = "Literature  [%]")

library(gridExtra)
grid.arrange(Starch, Protein, Oil, Fiber,Ash, nrow = 3)

