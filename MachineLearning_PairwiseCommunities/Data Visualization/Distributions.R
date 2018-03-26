# 18 Aug 17
# Demetrius DiMucci

# Figures 2B and 4B, 2C and 4C
# Distribution of growth responses for E.coli and for FBA

# E.coli responses
Ixns <- read.table('Interactions.txt',header=TRUE)
response <- log2(sort(Ixns[,5]))
pdf(file='EcoliDist.pdf')
par(mar=c(6,6,6,6))
plot(response,ylab=expression(paste('Log'[2],' Fold Change')),
     cex=2,cex.lab=2.5,ylim=c(-4,6.5),cex.axis=2,pch=1)

points(c(89,90),c(1,1),col='red',pch=16,cex=1)
points(rep(-4,9),cex=2,pch=16,col='darkgreen')
abline(h=1,lty=3,col='red',lwd=4)
legend('bottomright',legend=c('FC = 2','FC = 0'),pch=c(16,16,16),
       bg='white',col=c('red','darkgreen'),cex=1.4)
dev.off()

pdf(file='CostCorrelation.pdf')
par(mar=c(6,6,6,6))
R2 <- round(cor(Ixns[,7],Ixns[,5])^2)

plot(Ixns[,7],Ixns[,5],xlab='Biosynthetic Cost',ylab='Receiver Fold Change',
     cex.axis=2,cex=2,cex.lab=2.5)
legend('topright',legend=(expression(paste('R'^2,' .038'))),
       cex=1.4)
dev.off()

negEcoli <- vector()
count <- 0
for(i in 1:14){
  start <- 13*(i-1)+1
  end <- 13*i
  
  negEcoli[i] <- length(which(Ixns[start:end,5] <= 2))
  count <- count + 1
}

# Fig 4C
pdf(file='EcoliHist.pdf')
par(mar=c(5,5,5,5))
hist(negEcoli/13,col='grey',breaks=10,xlab='Weak Fraction',
     ylab=expression(paste(italic('E.coli '),'Strains')),main='',cex.axis=2,cex.lab=2.5)
dev.off()

# FBA responses
Outcomes <- read.table("RYTable.txt")
Monos <- which(Outcomes[,1]==Outcomes[,2])
Interactions <- Outcomes[-Monos,]
yields = sort(Interactions[,3])
pdf(file='FBADist.pdf')
par(mar=c(6,6,6,6))

neutral <- which(yields==0)
colors <- rep('black',9900)
colors[neutral]='red'
plot(yields,cex=.5,ylab='Relative Yield',cex.lab=2.5,cex.main=2,col=colors,cex.axis=2)
#abline(h=0,lwd=2,col='red',lty=2) 
legend('bottomright','Relative Yield = 0',pch=16,col='red',lwd=0,cex=1.4)
dev.off()



# This was not used in the text.
Jaccard_dist <- read.csv("JaccardPredictions.csv")[,2]

plot(Jaccard_dist,Interactions[,3])

pdf(file='JaccardFBACorrelation.pdf')
par(mar=c(6,6,6,6))
R2 <- round(cor(Interactions[,3],Jaccard_dist)^2,3)

plot(Jaccard_dist,Interactions[,3],xlab='Jaccard Distance',ylab='Relative Yield',
     cex.axis=2,cex=.7,cex.lab=2.5)
legend('topright',legend=(expression(paste('R'^2,' .012'))),
       cex=1.4)
dev.off()

negIxns <- vector()
accuracy <- vector()
count <- 0
for(i in 1:100){
  start <- 100*(i-1)+1
  end <- 100*i
  
  negIxns[i] <- length(which(OBS[start:end,3] < 0))
  count <- count + 1
}

accuracy <- vector()
for(i in 1:100){
  start <- 99*(i-1)+1
  end <- 99*i
  
  accuracy[i] <- length(which(exRf$y[start:end] == exRf$predicted[start:end]))
}

# Fig 2C in the text
pdf(file='FBADist_hist.pdf')
par(mar=c(5,5,5,5))
hist(negIxns/99,ylab='# Organisms',cex.axis=2,main='',
     cex.lab=2.5,cex=1,pch=16,xlab='Negative Fraction',col='grey')
dev.off()
####
# count <- 0
# corEcoli <- vector()
# for(i in 1:14){
#   start <- 13*(i-1)+1
#   end <- 13*i
#   
#   corEcoli[i] <- cor(Ixns[start:end,5],Ixns[start:end,8],method='spearman')
#   plot(Ixns[start:end,5],Ixns[start:end,8])
#   abline(0,1)
#   count <- count + 1
# }
