library(randomForest)
library(AUC)

load('EcoliRF')
exchanges <- read.table('ExchangeRxns.csv')
exchanges[,1] <- gsub('_p','.p',exchanges[,1])
# Biosynthetic costs
Ratio <- read.csv('Ratio_predict.csv')[,2]
Difference <- read.csv("Difference_predict.csv")[,2]
costs <- c(23,46,10,33,29,27,22,31,16,23,10,17,63,44)

responder <- vector()
for(i in 1:14){
  responder <- c(responder,rep(costs[i],13))
}
partner <- rep(costs,13)

pdf(file='EcoliROC.pdf')
par(mar=c(5,5,5,5))
plot(roc(mod1$votes[,2],mod1$y),col='black',cex.lab=2.5,cex.main=1.5,lwd=3,cex.axis=2)
#plot(roc(Ratio,mod1$y),col='darkgreen',add=TRUE)
plot(roc(Difference,mod1$y),col='blue',add=TRUE,lwd=3)
plot(roc(responder,mod1$y),col='red',add=TRUE,lwd=3)
plot(roc(partner,mod1$y),col='darkgreen',add=TRUE,lwd=3)

rfAuc = round((auc(roc(mod1$votes[,2],mod1$y))),3)
#rfRatio = round(auc(roc(Ratio,mod1$y)),3)
rfDif = round(auc(roc(Difference,mod1$y)),3)
rfResponder = round(auc(roc(responder,mod1$y)),3)
rfPartner = round(auc(roc(partner,mod1$y)),3)
legend('bottomright',
       legend=paste(c('Random Forest: ','Receiver Cost','Cost Difference','Giver Cost'),
                    c(rfAuc,rfResponder,rfDif,rfPartner)),
       col=c('black','red','blue','darkgreen'),pch=16,cex=1.4)
dev.off()

# # FBA ROC Curve using all features
# NoCV <- list.files('RF_fullData/',full.names=TRUE)
# load(NoCV[5])
Jaccard_dist <- read.csv("JaccardPredictions.csv")[,2]
# rf <- train.rf 
# 
# pdf(file='FBA_ROC.pdf')
# par(mar=c(5,5,5,5))
# plot(roc(train.rf$votes[,2],train.rf$y),cex.lab=1.5,cex.main=1.5)
# plot(roc(J_acc,Labels),col='red',add=TRUE)
# RFauc <- round(auc(roc(train.rf$votes[,2],train.rf$y)),3)
# JAuc = round((auc(roc(J_acc,Labels))),3)
# legend('bottomright',legend=paste(c('RF AUC: ','Jaccard AUC: '),c(RFauc,JAuc)),col=c('black','red'),pch=16)
# dev.off()

# FBA ROC Curve Exchanges only vs All

load('RF_exRxns')

pdf(file='FBA_ROC.pdf')
par(mar=c(5,5,5,5))
plot(roc(exRf$votes[,2],exRf$y),col='black',cex.lab=2.5,lwd=3,cex.axis=2)
exRFauc <- round(auc(roc(exRf$votes[,2],exRf$y)),3)
plot(roc(Jaccard_dist,exRf$y),col='red',add=TRUE,lwd=3)
Jauc <- round((auc(roc(Jaccard_dist,exRf$y))),3)
legend('bottomright',c(paste(c('Random Forest: '),c(exRFauc)),paste(c('Jaccard Distance: '),c(Jauc))),col=c('black','red'),pch=16,cex=1.4)
dev.off()

# Variable importance plots
rownames(mod1$importance) <- gsub('.2','.p',rownames(mod1$importance))
varImpPlot(mod1,type='1',main='',cex.lab=1)
rownames(exRf$importance)=exchanges[,1]
varImpPlot(exRf,type='1')

p=(mod1$importance[,3]/mod1$importanceSD[,3])
plot(sort(p),seq(1:28))

pdf(file='Ecoli_ranks.pdf')
varImpPlot(mod1,type='1',main='Auxotroph E.coli Variable Ranks')
dev.off()

pdf(file = 'GutRanks.pdf')
varImpPlot(exRf,type='1',main='')
dev.off()


####
coli=(mod1$importance[,3]/mod1$importanceSD[,3])
pdf(file='Ecoli.pdf',width=10,height=10)
par(mar=c(6,6,5,5))
plot(sort(coli),1:28, yaxt = "n",  ylab = "",
     xlab='',cex=2,cex.lab=3.5,cex.axis=2,
     font.lab=1, pch=16)
mtext('Mean Decrease Accuracy',side=1,line=4,cex = 4)
axis(2,seq(1:28),aminos[order(p)],las=2,cex.axis=2,font.axis=2)
for(i in 1:28){
  abline(h=i,lty=2,lwd=.5)
}
dev.off()

v = read.csv('Discovery_Table.csv')
metNames <- cbind(exchanges,'')
nameO <- rep(0,388)
for(i in 1:nrow(v)){
  found <- grep(v[i,1],metNames[,1])
  if(length(found )> 0){
    count <- count + 1
    print(found)
    nameO[found[1]] = v[i,2]
    nameO[found[2]] = paste0(v[i,2],'.p')
  }
}

exchanges[fba$ix[369:388],1] = c('D-Ribose','Galactose','L-Lysine','L-Glutamate.p','Betaine.p','L-Proline.p','GlcNac.p',
                               'Thiamin.p','dhptd','L-Glutamate','Thyminose.p','Thyminose','XAN.p',
                               'L-Arginine','L-Arabinose','XAN','D-Ribose.p','L-Lysine.p','dhptd.p','DAP.p')

rownames(exRf$importance)=nameO
pdf(file='FBA.pdf',width=12,height=10)
par(mar=c(6,13,5,5))
fba=(exRf$importance[,3]/exRf$importanceSD[,3])
fba = sort(fba,index.return=T)
plot(fba$x[369:388],1:20, yaxt = "n",  ylab = '',
     xlab='',cex=2,cex.lab=2,cex.axis=2,
     font.lab=2, pch=16)
mtext('Mean Decrease Accuracy',side=1,line=4,cex = 4)
axis(2,seq(1:20),exchanges[fba$ix[369:388],1],las=1,cex.axis=2,font.axis=2)
for(i in 1:20){
  abline(h=i,lty=2,lwd=.5)
}
dev.off()
