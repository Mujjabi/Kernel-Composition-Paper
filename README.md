---
title: "Dataset and Phenotypic Analysis Workflow - Kernel Composition Paper"
author: "Christopher Mujjabi"
date: "2023-04-28"
output: pdf_document
---
## Datasets Used in the Analysis  

### File "GrainQualityData_AmesPanel.csv" : 
Contains our raw data obtained from the NIR analysis of the inbred lines in the Ames diversity panel 

### File "Wis-Rawdata.csv": 
Contains part of the raw NIR data from Renk 2021 supplimental tables, for a set of 275 common genotypes between Renk 2021 study and our study. This dataset contains the 275 genotypes tested in 5 environments with 2 replicatons per environment.

### File "TableValue.csv": 
Contains the compositional trait means and estimated Blups for each heterotic maize group obtained in our study (Ames_means & Ames_Blups) and the trait Blups of the similar heterotic groups reported in the Renk 2021 paper (Wisc_Blups). The values in the Wisc_Blups were obtained from Table 2 (Summary of best linear unbiased prediction (BLUP)-adjusted mean differences for each compositional trait compared among maize types) in the Renk 2021 paper. TableValue file was used to conduct the correlations between the obtained Blups in the 2 studies. 

## Phenotypic Analysis Workflow 

### File "Final Analysis_AmesPaper.R":
Contains the R Code used to conduct the phenotypic analysis and comparison of data between replicared (Renk 2021) and our unreplicated experiment. 

### R Code and Workflow Explained
This is the final workflow for the phenotypic analysis of kernel compositional traits.

##### R Packages and Libraries Used
Here are the packages and libraries used in the analysis. 

```{r}
library(sf)
library(dplyr)
library(ggplot2)
library(tmap)
```










First, I conducted an analysis of variance between heterotic groups since this wasnt possible between genototypes due to lack of replications in our experiment.I used a simple linear model which assumes fixed factors, since we dont have random factors. Heterotic groups are  fixed groups thus considered fixed factors in the mode. 

```{r setup, include=FALSE, echo = TRUE, warning = F, message = F}

```


