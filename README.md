# Dataset and Phenotypic Analysis Workflow - Kernel Composition Paper"
## Author: Christopher Mujjabi
## Date: 2023-04-28

### Datasets Used in the Analysis  

#### File "GrainQualityData_AmesPanel.csv" : 
Contains our raw data obtained from the NIR analysis of the inbred lines in the Ames diversity (NCRPIS) panel.

#### File "Wis-Rawdata.csv": 
Contains part of the raw NIR data from Renk 2021 supplimental tables, for a set of 275 common genotypes between Renk 2021 study and our study. This dataset contains the 275 genotypes tested in 5 environments with 2 replicatons per environment.

#### File "TableValue.csv": 
Contains the compositional trait means and estimated Blups for each heterotic maize group obtained in our study (Ames_means & Ames_Blups) and the trait Blups of the similar heterotic groups reported in the Renk 2021 paper (Wisc_Blups). The values in the Wisc_Blups were obtained from Table 2 (Summary of best linear unbiased prediction (BLUP)-adjusted mean differences for each compositional trait compared among maize types) in the Renk 2021 paper. TableValue file was used to conduct the correlations between the obtained Blups in the 2 studies. 

#### File "Final Analysis_AmesPaper.R":
Contains the R Code used to conduct the phenotypic analysis and comparison of data between replicared (Renk 2021) and our unreplicated experiment.

### Phenotypic Analysis Workflow and R Code Explained
This is the final workflow for the phenotypic analysis of kernel compositional traits.

#### 1. R Packages and Libraries Used
Here are the packages and libraries used in the analysis. 

```{r}
library(agricolae)
library(lme4)
library(jtools)
library(lmerTest)
library(car)
library(dplyr)
library(Hmisc)
library("ggpubr")
library(gridExtra)
library(xlsx)
library(readr)
```
#### 2. Uploading the dataset
First, I uploaded the dataset we obtained from scanning a sub-set of 954 inbred lines in the Ames diversity panel.  This file is called "GrainQualityData_Amespanel.csv". Then I changed each column in the dataset to either "factor" if categorical or "numeric" if numeric variable. 

```{r}
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
```
#### 3. Data Cleaning

Before conducting any statistical procedures, I first cleaned the data to detect and remove outliers and any data points that seemed off. First, I plotted boxplots to see the distribution of each trait and detect any outliers within each trait. 
```{r}
boxplot(Compiled[4:9], plot=TRUE)$out
```
Then, I used the interquartile range rule (IQR method) to identify and eliminate outliers for each trait using the code below. 

```{r}
outliers <- function(x) {
  
  Q1 <- quantile(x, probs=.25) 
  Q3 <- quantile(x, probs=.75) 
  iqr = Q3-Q1           
  
  upper_limit = Q3 + (iqr*1.5) 
  lower_limit = Q1 - (iqr*1.5)  
  
  x > upper_limit | x < lower_limit  
}

remove_outliers <- function(Compiled, cols = names()) {  
  for (col in cols) {
    Compiled <- Compiled[!outliers(Compiled[[col]]),] 
  }
  Compiled
}

Compiled2 <- remove_outliers(Compiled, c(4:9))

outliers <- boxplot(Compiled[4:9], plot=TRUE)$out  
x<-Compiled
Compiled <- x[-which(x$Starch %in% outliers),] 
Compiled <- Compiled[-which(Compiled$Protein %in% outliers),]
Compiled <- Compiled[-which(Compiled$Oil %in% outliers),]
Compiled <- Compiled[-which(Compiled$Fib %in% outliers),]
Compiled <- Compiled[-which(Compiled$Ash %in% outliers),]
Compiled <- Compiled[-which(Compiled$Density %in% outliers),]
```
A clean dataset without outliers was created and named "Compiled2". Here is a new box plot showing no outliers within each kernel composition trait. 
```{r}
boxplot(Compiled2[4:9], col=c("red", "yellow", "green","#999999", "#E69F00", "#56B4E9"), plot=TRUE)$out 
```
#### 4. Analysis Of Variance (ANOVA)
First, I conducted an analysis of variance (ANOVA) between heterotic groups since this wasnt possible between genotypes due to lack of replications in our experiment.I used a simple linear model which assumes fixed factors, since we dont have random factors. Heterotic groups are  fixed groups thus considered fixed factors in the mode. 

```{r setup, include=FALSE, echo = TRUE, warning = F, message = F}
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
```
#### 5. Post-Hoc Test (Tukeys-HSD Test)
I used a Tukeys-HSD Test to conduct a pairwise mean comparisons for each kernel composition trait, within each heterotic group. This can be done in 2 ways as showed in the codes below. 
```{r setup, include=FALSE, echo = TRUE, warning = F, message = F}
library(agricolae)
Protein <-HSD.test(modelprot,"Group", group=TRUE)
print(Protein)

Starch  <-HSD.test(modelstr,"Group", group=TRUE)
print(Starch)

Oil <-HSD.test(modeloil,"Group", group=TRUE)
print(Oil)

Ash <-HSD.test(modelash,"Group", group=TRUE)
print(Ash)

Fiber  <-HSD.test(modelfiber,"Group", group=TRUE)
print(Fiber)

Density <-HSD.test(modeldensity,"Group", group=TRUE)
print(Density)
```
OR 
```{r setup, include=FALSE, echo = TRUE, warning = F, message = F}
df<-df.residual(modelprot)
MSerror<-deviance(modelprot)/df
with(Compiled2,HSD.test(Protein,Group,df,MSerror, group=TRUE,console=TRUE,
                          main="Protein Composition for Different Maize Groups "))
```
#### 6. Estimating BLUPs for Each Trait 

I used the lmer function to fit a mixed linear model, using Heterotic groups as random variables to be able to extract Best Linear Unbiased predictors (BLUPs)-adjusted trait means for each heterotic group. Blups were extracted using the "coef" function. The function calculates estimated variances between random-effects terms in a mixed-effects model 

```{r setup, include=FALSE, echo = TRUE, warning = F, message = F}
mod1 <- lmer(Protein ~ (1|Group), REML = TRUE, data= Compiled)
mod2 <- lmer(Starch ~ (1|Group), REML = TRUE, data= Compiled)
mod3 <- lmer(Oil ~ (1|Group), REML = TRUE, data= Compiled)
mod4 <- lmer(Ash ~ (1|Group), REML = TRUE, data= Compiled)
mod5 <- lmer(Fib ~ (1|Group), REML = TRUE, data= Compiled)
mod6 <- lmer(Density ~ (1|Group), REML = TRUE, data= Compiled)

varComp<-as.data.frame(VarCorr(mod1,comp="vcov")) 
blupProt= coef(mod1)$Group

varComp<-as.data.frame(VarCorr(mod2,comp="vcov"))
blupStr = coef(mod2)$Group

varComp<-as.data.frame(VarCorr(mod3,comp="vcov"))
blupOil = coef(mod3)$Group

varComp<-as.data.frame(VarCorr(mod4,comp="vcov"))
blupash = coef(mod4)$Group

varComp<-as.data.frame(VarCorr(mod5,comp="vcov"))
blupfiber = coef(mod5)$Group

varComp<-as.data.frame(VarCorr(mod6,comp="vcov"))
blupdensity = coef(mod6)$Group

Ames_Blups <- cbind(blupProt,blupStr,blupOil,blupfiber,blupdensity,blupash)

write.csv(Ames_Blups, "Ames_Blups.csv")
```
The obtained BLUPs above were combined with the BLUPs reported in the Renk et al, 2021 paper (Table 2) and formed a new dataset called "TableValues". The columns in this dataset are: Heterotic Group, Trait, Ames_Blups and Wisc_Blups for each trait for each trait and was used in the next step to perform correlation analysis. 

##### 6.1 Test for normality Assumptioms
I used the qq-plots to test for any violations of the model assumptions (independence of residuals, residual normality, and homoscedastic variance of residuals) using the line of code below.

```{r setup, include=FALSE, echo = TRUE, warning = F, message = F}
par(mfrow=c(2,3)

qqPlot(residuals(mod1), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Starch Content")

qqPlot(residuals(mod2), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Protein Content") 

qqPlot(residuals(mod3), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Oil Content") 

qqPlot(residuals(mod4), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Ash Content") 

qqPlot(residuals(mod5), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Fiber Content") 

qqPlot(residuals(mod6), pch=19, col="dark blue", col.lines="red", xlab="Predicted quantiles", ylab="Observed quantiles", main = "Kernel Density") 
```
#### 7. Correlation Analysis

##### 7.1 Spearman Ranking Correlation 
First, I conducted a spearman ranking correlation to compare how the heterotic groups performed in the two studies. Strong correlations would validate the data obtained from our unreplicated experiment. 

```{r setup, include=FALSE, echo = TRUE, warning = F, message = F}
Data1  <- read.csv("TableValue.csv", header = TRUE)

library(dplyr) 
Corr1 <- Data1%>% 
  group_by (Trait) %>% 
  summarise(cor=cor(Ames_Blups,Wisc_Blups, method = "spearman"))

write.csv(Corr1, "spearmanCorr.csv")
```
##### 7.2 Pearson Correlations 

In addition, we conducted a Pearson correlation between the kernel composition trait values obtained from the Renk et al, 2021 replicated experiment and values obtained from our unreplicated experiment, for the 275 common genotypes between the 2 studies. The raw data for the Renk et al, 2021 experiment was downloaded from the supplimentary tables. The common genotypes were identified, and the rest of the genotyopes were deleted from the dataset. This file was names as "Wis-Rawdata.csv". The dataset has 275 genotypes, tested in 5 environments with 2 replicates in each environment. The dataset was uploaded and trait means for each genotype were estimated in the code below. 

```{r setup, include=FALSE, echo = TRUE, warning = F, message = F}
Wis <- read.csv("Wis-Rawdata.csv", header = TRUE)  
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
  summarise(across(Prot_W:Star_W, mean)) 
write.csv (WisMeans, "WisMeans.csv") 
```
The estimated trait means from the Renk et al, 2021 study were Combined with values obtained from our unreplicated experiment using the "merge" function, to form a single dataset called "WisAmesCombined" as showed in the code below.

```{r setup, include=FALSE, echo = TRUE, warning = F, message = F}
WisAmesCombined <- merge(Compiled, WisMeans, by = 'Accession') #merge function joins the 2 datasets by the common accession 
write.csv (WisAmesCombined, "Combined2.csv") 
```
Using the generated dataset (WisAmesCombined),a Pearson Correlation between our observed trait values from an unreplicated experiment and the Renk et al, 2021 trait meav values from a replicated experiment, for the 275 common genotypes, was conducted using the code below. 

```{r setup, include=FALSE, echo = TRUE, warning = F, message = F}
library(Hmisc) 
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
write.csv (Corr3, "PearsonCorr.csv")
```
The Generated file Corr3 shows the correlation between the trait values for each genotype. 

##### 7.3 Correlation Scatter Plots

The obtained pearson correlations were plotted on scatter plots using the ggscatter function as shown below. 

```{r setup, include=FALSE, echo = TRUE, warning = F, message = F}
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
```



