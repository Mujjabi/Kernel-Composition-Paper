
# Dataset and Phenotypic Analysis Workflow - Kernel Composition Paper"
## Author: Christopher Mujjabi
## Date: 2023-04-28

### Datasets Used in the Analysis  

#### File "GrainQualityData_AmesPanel.csv" : 
Contains our raw data obtained from the NIR analysis of the inbred lines in the Ames diversity panel 

#### File "Wis-Rawdata.csv": 
Contains part of the raw NIR data from Renk 2021 supplimental tables, for a set of 275 common genotypes between Renk 2021 study and our study. This dataset contains the 275 genotypes tested in 5 environments with 2 replicatons per environment.

#### File "TableValue.csv": 
Contains the compositional trait means and estimated Blups for each heterotic maize group obtained in our study (Ames_means & Ames_Blups) and the trait Blups of the similar heterotic groups reported in the Renk 2021 paper (Wisc_Blups). The values in the Wisc_Blups were obtained from Table 2 (Summary of best linear unbiased prediction (BLUP)-adjusted mean differences for each compositional trait compared among maize types) in the Renk 2021 paper. TableValue file was used to conduct the correlations between the obtained Blups in the 2 studies. 

### Phenotypic Analysis Workflow 

#### File "Final Analysis_AmesPaper.R":
Contains the R Code used to conduct the phenotypic analysis and comparison of data between replicared (Renk 2021) and our unreplicated experiment. 

#### R Code and Workflow Explained
This is the final workflow for the phenotypic analysis of kernel compositional traits.

##### R Packages and Libraries Used
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
```
##### Uploading the dataset
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
##### Data Cleaning

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
Here is a new box plot showing no outliers within each kernel composition trait. 
```{r}
boxplot(Compiled2[4:9], plot=TRUE)$out 
```
##### Analysis Of Variance (ANOVA)
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
##### Post-Hoc Test (Tukeys-HSD Test)
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




