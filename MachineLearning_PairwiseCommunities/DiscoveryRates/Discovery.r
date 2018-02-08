# 8 Jan 17
# Demetrius DiMucci
library(forestFloor)
library(randomForest)
Data <- read.table('ConcatenatedVectors.txt')
Outcomes <- read.table("RYTable.txt")
Monos <- which(Outcomes[,1]==Outcomes[,2])

Elements <- (read.table('Curated_Element_Names.txt'))
Elements <- (Elements[,1])
Elements <- as.vector(Elements)
JointElements <- c(Elements,Elements)
JointElements[2084:4166] <- paste0(JointElements[2084:4166],'.p')

data <- Data[-Monos,]
use <- grep('EX',JointElements)
exData <- data[,use]

load('RF_exRxns')
ff = forestFloor(exRf,exData,binary_reg = T,calc_np=F,bootstrapFC = F)
X <- ff$FCmatrix

## Sample 1 the first mechanism was encountered at combined rank 1 and second at
## combined rank 2. The top mechanism was variable 89
colors <- rep('grey',17)
colors[7] <- 'red'
begin = 165;
end = 180;
index1 <- c(81,begin:end)
begin1 = begin + 194;
end1 <- end + 194
index2 <- c(81+194,begin1:end1)
pdf(file='FirstHalf.pdf')
barplot(X[2,index1],col=colors,names='',yaxt='n')
dev.off()
pdf(file='SecondHalf.pdf')
barplot(X[2,index2],col=colors,names='',yaxt='n')
dev.off()
Y <- X[2,index1]+X[2,index2]
j <- sort(Y,index.return=T)
pdf(file='Combined.pdf')
barplot(Y,col=colors,names='',yaxt='n')
dev.off()
pdf(file='Combined_sorted.pdf')
barplot(sort(Y),col=colors[j$ix],names='',yaxt='n')
dev.off()
#####

Positions <- read.table('Positions.csv')
negs <- which(Positions[,12] == -1)
pos <- which(Positions[,12] == 1)
realPos <- which(Outcomes[-Monos,3]>0)

# Columns of Positions explained

# FirstNeg (col 1) - rank where the first negative mechanism was encountered
# SecondNeg (col 2) - rank where the second negative mechanism was encountered, if it exists
# FirstPos (col 3) - rank where the first positive mechanism was encountered
# SecondPos (col 4) - rank where the second positive mechanism was encountered, if it exists
# TotalMechs (col 5) - total number of metabolites both models have a transport reaction for, really is potential mechanisms
# TotalNeg (col 6) - total number of metabolites that were subject to competition in the simulation
# TotalPos (col 7) - total number of metabolites subject to cross feeding in the simulation
# FalseMech1 (col 8) - total number of potential metabolites (col 6) encountered before encountering a real one
# FalseMech2 (col 9) - total number of potential metabolites (col 6) encountered before encountering a second real one
# FalseMech3 (col 10) - total number of potential metabolites (col 6) encountered before encountering a real one, starting from position 194
# FalseMech4 (col 11) - total number of potential metabolites (col 6) encountered before encountering a second real one, starting from position 194


NEGS <- Positions[negs,]
POS <- Positions[pos,]
POS1 <- Positions[realPos,]
# calculate the expected position of the first encounter for each value of the total number of 
# real mechanisms
expected <- vector()
#expected2 <- vector()
for(i in 1:max(Positions[,5])){
  b = 194 - i
  expected[i] <- (b + i + 1)/(i + 1)
}

# Uncommnet this to see expected encounter rates if we looked only at 
# Potential mechanisms - those where both models have the transporter
expectedSharedOnly <- vector()
for(i in 1:nrow(NEGS)){
  real = NEGS[i,6]
  b = NEGS[i,5] - real
  expectedSharedOnly[i] <- (b + real + 1)/(real + 1)
}

expectedSharedOnlyPos <- vector()
for(i in 1:nrow(POS1)){
  real = POS1[i,7]
  b = POS1[i,5] - real
  expectedSharedOnlyPos[i] <- (b + real + 1)/(real + 1)
}

# Determine the expected distribtion of ranks for this data set
# encounterExpected is the expected rank of the first competitive metabolite to be found
encounterExpected <- vector()
for(i in 1:nrow(NEGS)){
  encounterExpected[i] <- expected[NEGS[i,6]]
}

posExpected <- vector()
for(i in 1:nrow(POS1)){
  posExpected[i] <- expected[POS1[i,7]]
}

foundPos <- vector()
for(i in 1:420){
  foundPos[i] <- rank(j)[POS1[i,3]]
}

pdf(file='ExpectedRanks.pdf')
hist(encounterExpected,breaks=200,xlab='Rank',main = 'Expected Rank of first discovered metabolite')
dev.off()

pdf(file='FirstMechanism.pdf')
# Remove 2 outliers for purposes of visibility
hist(NEGS[-c(2570,4197),1],breaks=200,main='Rank of first discovered metabolite',xlab='Rank')
dev.off()

# Visualize where 1 and 2 mechanism samples are found
ones <- which(NEGS[,6] == 1)
twos <- which(NEGS[,6] == 2)

pdf(file='OneMechanism.pdf')
# Remove 2 outliers for purposes of visibility
hist(NEGS[ones,1],breaks=200,main='Rank of first discovered metabolite',xlab='Rank')
mtext('One Mechanism')
dev.off()

pdf(file='TwoMechanisms.pdf')
# Remove 2 outliers for purposes of visibility
hist(NEGS[twos,1],breaks=200,main='Rank of first discovered metabolite',xlab='Rank')
mtext('Two Mechanisms')
dev.off()

# Density Plots
hist(NEGS[,1],freq = FALSE,main='',xlab='',ylab = '',axes=T,col='grey',ylim=c(0,.35),breaks=200)
lines(density(NEGS[,1]))

pdf("ExperimentsNeeded.pdf")
par(mar=c(5,5,5,5))
plot(density(NEGS[,1]),xlab='Experiments Before 1st Discovery',main=''
     ,cex.lab=2,cex.axis=2,ylab='Probability')
N <- density(NEGS[,1])
x <- N$x
y <- N$y
polygon(c(x,-1),c(y,0),col='black')
points(density(encounterExpected),type='l',col='red',lwd=3)
legend('topright',c('Using Feature Contributions','Random'),pch=c(16,16),
       col=c('black','red'))
dev.off()


pdf('OneMechanism.pdf')
par(mar=c(5,5,5,5))
ones <- which(NEGS[,6]==1)
plot(density(NEGS[ones,1]),xlab='Experiments Before 1st Discovery',main='',
     cex.lab=2,cex.axis=2,ylab='Probability')
N <- density(NEGS[ones,1])
x <- N$x
y <- N$y
polygon(c(x,-2),c(y,-0.0),col='black')
segments(0,1/194,185,1/194,lwd=3,col='red')
legend('topright',c('Using Feature Contributions','Random'),pch=c(16,16),
       col=c('black','red'))
dev.off()


####
pdf("ExperimentsNeeded_PositveCases.pdf")
par(mar=c(5,5,5,5))
plot(density(foundPos),xlab='Experiments Before 1st Discovery',main=''
     ,cex.lab=2,cex.axis=2,ylim=c(0,.03),ylab='Probability')
N <- density(foundPos)
x <- N$x
y <- N$y
polygon(c(x,-1),c(y,0),col='black')
points(density(posExpected),type='l',col='red',lwd=3)
legend('topright',c('Using Feature Contributions','Random'),pch=c(16,16),
       col=c('black','red'))
dev.off()

# Make a table to identify which metabolites are used as strong predictors
library(randomForest)
library(forestFloor)
load('RF_exRxns')

Data <- read.table('ConcatenatedVectors.txt')
Outcomes <- read.table('RYTable.txt')
labels <- Outcomes[,3]
labels[which(labels < 0)] = -1
labels[which(labels >= 0)] = 1
Monos <- which(Outcomes[,1]==Outcomes[,2]) # Monos are those samples that were monoculture


Elements <- (read.table('Curated_Element_Names.txt'))
Elements <- (Elements[,1])
Elements <- as.vector(Elements)
JointElements <- c(Elements,Elements)
JointElements[2084:4166] <- paste0(JointElements[2084:4166],'.1')

data <- cbind(Data,labels)
data <- data[-Monos,]

use <- grep('EX',JointElements)
exData <- data[,use]
Labels <- as.factor(data[,ncol(data)])
# if you want to train a new random forest uncomment this
# exRF_new <- randomForest(Labels ~.,data=exData)

# Forest floor

ff = forestFloor(exRf,exData,binary_reg = T,calc_np=F,bootstrapFC = F)

# Store the feature contribution matrix
X <- ff$FCmatrix

# colnames(ff$X)=JointElements[use]
# Col = fcol(ff,cols=20)
# pdf(file='GOV_FBA.pdf')
# plot(ff,col=Col,plot_GOF = T,limitY = T,speedup_GOF = T,plot_seq = 1:28)
# dev.off()
# The positions table tells us the index of the first mechanism identified
# We just rank each sample again and look at which sample is rank X
Interactions <- Outcomes[-Monos,]

metabolites <- vector()
respondingSpecies <- vector()
partnerSpecies <- vector()
compounds <- JointElements[use[1:194]]

Xnegs <- X[negs,]
Xixns <- Interactions[negs,]
for(i in 1:nrow(NEGS)){
  tempRank <- sort((Xnegs[i,1:194] + Xnegs[i,195:388]),index.return=T)
  firstMech <- NEGS[i,1]
  metIndex <- tempRank$ix[firstMech]
  metabolites[i] <- compounds[metIndex]
  respondingSpecies[i] <- Xixns[i,1]
  partnerSpecies[i] <- Xixns[i,2]
}


MetStats <- matrix(0,nrow=length(unique(metabolites)),ncol=7)
rownames(MetStats) = unique(metabolites)
for(i in 1:nrow(MetStats)){
  met <- rownames(MetStats)[i]
  Ndex <- which(metabolites == met)
  MetStats[i,1] <- length(Ndex)
  MetStats[i,2] <- length(unique(respondingSpecies[Ndex]))
  MetStats[i,3] <- length(unique(partnerSpecies[Ndex]))
  MetStats[i,4] <- median(NEGS[Ndex,1])
  MetStats[i,5] <- median(NEGS[Ndex,8])
  MetStats[i,6] <- median(NEGS[Ndex,5])
  MetStats[i,7] <- median(NEGS[Ndex,6])
}

MetStats <- MetStats[rev(order(MetStats[,1])),]
colnames(MetStats) <- c('Times Found First','Unique Responders','Unique Partners','Median Rank','Median False Hits','Median Potential Mechs','Median Mechs')

rownames(MetStats) = gsub('EX_','',rownames(MetStats))
rownames(MetStats) = gsub('_e0','',rownames(MetStats))
rownames(MetStats) = gsub('_c0','',rownames(MetStats))

write.csv(MetStats,file='Discovery_Table.csv')


# Make table 1 - summary statistics for the two classes
SummaryTable <- matrix(0,ncol=6,nrow=2)
colnames(SummaryTable) = c('Samples','Median Potential Interaction Mets','Max Potential Mets','Min Potential Mets','Competitve Mets','CrossFeed Mets')
rownames(SummaryTable) <- c('Neg RY','Non-neg RY')
SummaryTable[1,1] = nrow(NEGS)
SummaryTable[1,2] = median(NEGS[,5])
SummaryTable[1,3] = max(NEGS[,5])
SummaryTable[1,4] = min(NEGS[,5])
SummaryTable[1,5] = median(NEGS[,6])
SummaryTable[1,6] = median(NEGS[,7])

SummaryTable[2,1] = nrow(POS)
SummaryTable[2,2] = median(POS[,5])
SummaryTable[2,3] = max(POS[,5])
SummaryTable[2,4] = min(POS[,5])
SummaryTable[2,5] = median(POS[,6])
SummaryTable[2,6] = median(POS[,7])

write.csv(SummaryTable,file='SummaryTable.csv')
  

#
nums <- sort(unique(NEGS[,6]))

FCbetter <- vector()
RankBetter <- vector()
for(i in 1:length(nums)){
  mechs <- which(NEGS[,6]==nums[i])
  
  FCbetter[i] <- length(which(H[mechs] < Y[mechs]))/length(mechs)
  RankBetter[i] <- mean(H[mechs]-Y[mechs])
}

nums <- sort(unique(realPos[,6]))
expectedTonly <- vector()
for(i in 1:nrow(realPos)){
  real = realPos[i,6]
  b = realPos[i,5] - real
  expectedTonly[i] <- (b + real + 1)/(real + 1)
}

FCbetter <- vector()
RankBetter <- vector()
for(i in 1:length(nums)){
  mechs <- which(realPos[,6]==nums[i])
  
  FCbetter[i] <- length(which(H[mechs] < Y[mechs]))/length(mechs)
  RankBetter[i] <- median(H[mechs]-Y[mechs])
}


#

mechFound <- vector()
guessed <- vector()
totMech <- vector()
mechs <- sort(unique(NEGS[,6]))

for(i in 1:length(unique(NEGS[,6]))){
  h <- which(NEGS[,6]==mechs[i])
  mechFound[i] <- median(NEGS[h,1])
  guessed[i] <- median(encounterExpected[h])
  totMech[i] <- length(h)
}


h <- 0
total=50000
for(i in 1:total){
  h <- h + sort((sample(1:194,2)))[2]
}
h/total
