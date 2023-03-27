#Libraries
library(Hmisc) #Describe Function
library(psych) #Multiple Functions for Statistics and Multivariate Analysis
library(GGally) #ggpairs Function
library(ggplot2) #ggplot2 Functions
library(vioplot) #Violin Plot Function
library(corrplot) #Plot Correlations
library(REdaS) #Bartlett's Test of Sphericity
library(psych) #PCA/FA functions
library(factoextra) #PCA Visualizations
library("FactoMineR") #PCA functions
library(ade4) #PCA Visualizations
##############################################################################################

#Set Working Directory
setwd('C:/Users/Samuel/Downloads/Advaced Data Analysis/Final Project')

#Read in Datasets
RawData <- read.csv(file="OnlineNewsPopularity.csv", header=TRUE, sep=",")

#Check Sample Size and Number of Variables
dim(RawData)
#39644-Sample Size and 61 variables

#Show for first 6 rows of data
head(RawData)

#Names of the data
names(RawData)

#Check for Missing Values (i.e. NAs)
#For All Variables
sum(is.na(RawData))
#0 total missing values (0 cells with missing data)

describe(RawData)

#Create new subsets of data 
Data <- RawData[,c(5:7,12,19:22,24:31,42,43,45:51,53,55,61)]

library(psych)
describe(Data)

#Create Initial Linear Regression Model with Enter Method
model1 <- lm(shares ~ ., data=Data)
model1

library(car)
#Check VIF
vif(model1)

# PCA_Plot functions
PCA_Plot = function(pcaData)
{
  library(ggplot2)
  
  theta = seq(0,2*pi,length.out = 100)
  circle = data.frame(x = cos(theta), y = sin(theta))
  p = ggplot(circle,aes(x,y)) + geom_path()
  
  loadings = data.frame(pcaData$rotation, .names = row.names(pcaData$rotation))
  p + geom_text(data=loadings, mapping=aes(x = PC1, y = PC2, label = .names, colour = .names, fontface="bold")) +
    coord_fixed(ratio=1) + labs(x = "PC1", y = "PC2")
}

PCA_Plot_Secondary = function(pcaData)
{
  library(ggplot2)
  
  theta = seq(0,2*pi,length.out = 100)
  circle = data.frame(x = cos(theta), y = sin(theta))
  p = ggplot(circle,aes(x,y)) + geom_path()
  
  loadings = data.frame(pcaData$rotation, .names = row.names(pcaData$rotation))
  p + geom_text(data=loadings, mapping=aes(x = PC3, y = PC4, label = .names, colour = .names, fontface="bold")) +
    coord_fixed(ratio=1) + labs(x = "PC3", y = "PC4")
}

PCA_Plot_Psyc = function(pcaData)
{
  library(ggplot2)
  
  theta = seq(0,2*pi,length.out = 100)
  circle = data.frame(x = cos(theta), y = sin(theta))
  p = ggplot(circle,aes(x,y)) + geom_path()
  
  loadings = as.data.frame(unclass(pcaData$loadings))
  s = rep(0, ncol(loadings))
  for (i in 1:ncol(loadings))
  {
    s[i] = 0
    for (j in 1:nrow(loadings))
      s[i] = s[i] + loadings[j, i]^2
    s[i] = sqrt(s[i])
  }
  
  for (i in 1:ncol(loadings))
    loadings[, i] = loadings[, i] / s[i]
  
  loadings$.names = row.names(loadings)
  
  p + geom_text(data=loadings, mapping=aes(x = PC1, y = PC2, label = .names, colour = .names, fontface="bold")) +
    coord_fixed(ratio=1) + labs(x = "PC1", y = "PC2")
}

PCA_Plot_Psyc_Secondary = function(pcaData)
{
  library(ggplot2)
  
  theta = seq(0,2*pi,length.out = 100)
  circle = data.frame(x = cos(theta), y = sin(theta))
  p = ggplot(circle,aes(x,y)) + geom_path()
  
  loadings = as.data.frame(unclass(pcaData$loadings))
  s = rep(0, ncol(loadings))
  for (i in 1:ncol(loadings))
  {
    s[i] = 0
    for (j in 1:nrow(loadings))
      s[i] = s[i] + loadings[j, i]^2
    s[i] = sqrt(s[i])
  }
  
  for (i in 1:ncol(loadings))
    loadings[, i] = loadings[, i] / s[i]
  
  loadings$.names = row.names(loadings)
  
  print(loadings)
  p + geom_text(data=loadings, mapping=aes(x = PC3, y = PC4, label = .names, colour = .names, fontface="bold")) +
    coord_fixed(ratio=1) + labs(x = "PC3", y = "PC4")
}

#Test KMO Sampling Adequacy
library(psych)
KMO(Data)
#Overall MSA =  0.6
#> 0.7 so good size

#Test Bartletts Test of Sphericity
library(REdaS)
bart_spher(Data)
#p-value < 2.22e-16 (Very Small Number)

#Test for Reliability Analysis using Cronbachs Alpha
library(psych)
alpha(Data,check.keys=TRUE)
#raw_alpha = 0.34

#scaling
scale(Data, center = TRUE, scale = TRUE)

#Create PCA
p = prcomp(Data, center=T, scale=T)
p

#Check Scree Plot
plot(p)
abline(1, 0)

#Check PCA Summary Information
summary(p)
print(p)

#Check PCA visualizations
plot(p) #Scree Plot
PCA_Plot(p) #PCA_plot1
PCA_Plot_Secondary(p) #PCA_Plot2
#biplot(p) #Biplot

#Calculating the Varimax Rotation Loadings manually
rawLoadings = p$rotation %*% diag(p$sdev, nrow(p$rotation), nrow(p$rotation))
print(rawLoadings)
v = varimax(rawLoadings)

#Options available under varimax function
ls(v)
v

# The Psych package has a wonderful PCA function that allows many more options
# including build-in factor rotation, specifiying a number of factors to include 
# and automatic "score" generation

#Best Way to Conduct PCA Analysis
p2 = psych::principal(Data, rotate="varimax", covar=FALSE, nfactors=8, scores=TRUE)
p2
print(p2$loadings, cutoff=.4, sort=T)

#PCAs Other Available Information
ls(p2)
p2$values
table(p2$values>1)
p2$communality
p2$rot.mat

#Calculating scores
scores <- p2$scores
scores_1 <- scores[,1]
round(cor(scores),2)
cor(scores)
summary(scores)
summary(scores_1)
scores_2 <- scores[,2]
summary(scores_2)
scores_3 <- scores[,3]
summary(scores_3)
scores_4 <- scores[,4]
summary(scores_4)
scores_5 <- scores[,5]
summary(scores_5)
scores_6 <- scores[,6]
summary(scores_6)
scores_7 <- scores[,7]
summary(scores_7)
scores_8 <- scores[,8]
summary(scores_8)
#scores_9 <- scores[,9]

