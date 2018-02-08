# 17 Jan 2018
# Demetrius DiMucci
# Create Barplots For Fig 5

library(forestFloor)
library(randomForest)
library(gridExtra)
library(grid)
library(gtable)

load('EcoliRF')
Ixns <- read.table("Interactions.txt",header=TRUE)
data <- read.table('Observations.txt',header=TRUE)
Y <- data[,-29]

# Define Column Names
aminos1 <- c('C','F','G','H','I','K','L','M','P','R','S','T','W','Y')
aminos2 <- paste0(aminos1,'.p')
aminos <- c(aminos1,aminos2)


ff = forestFloor(mod1,Y,binary_reg = T,calc_np=T,bootstrapFC = F)
X <- ff$FCmatrix
colnames(X) <- aminos

# Find which rows of X correspond to Proline with Serine
target <- which(Ixns[,3] == 'P' & Ixns[,4] == 'S')
pdf(file='Proline_Serine.pdf')
par(mar=c(5,5,5,5))
barplot(X[target,],names=aminos,ylim=range(X[target,]),las=3,ylab='Feature Contribution',
        cex.names= 1.25,cex.lab=2.5,cex.axis = 2)
dev.off()
# Find which rows of X correspond to Methionine with Cysteine
target <- which(Ixns[,3] == 'M' & Ixns[,4] == 'C')
pdf(file='Methionine_Cysteine.pdf')
par(mar=c(5,5,5,5))
barplot(X[target,],names=aminos,ylim=range(X[target,]),las=3,ylab='Feature Contribution',
        cex.names= 1.25,cex.lab=2.5,cex.axis = 2)
dev.off()
