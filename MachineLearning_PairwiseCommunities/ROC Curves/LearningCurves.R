library(Hmisc)

# Plot the E. coli learning curve
fractions_E <- c(.95,.9,.8,.7,.6,.5,.4,.3,.2,.1,.05)
medianRFs <- read.csv('MedianRF_ecoli.csv')[,2]
sdRFs <- read.csv('sdRFs_ecoli.csv')[,2]

pdf(file='EcoliCurve.pdf')
par(mar=c(5,5,5,5))
plot(fractions_E,medianRFs,ylim=c(.45,1),xlab='Training Set Fraction',type='l',
     ylab='Balanced Accuracy',cex.lab=2.5,cex.main=1.5,lwd=3,cex.axis=2)
minor.tick(nx=4,ny=2,tick.ratio=.5)
#points(fractions,medianRFs+sdRFs,col='red',cex=.5,type='l')
#points(fractions,medianRFs-sdRFs,col='red',cex=.5,type='l')
dev.off()

# Plot the FBA learning Curve
# Set the color palette
col_set=colorRampPalette(c("blue","red"))(90)

matrices <- list.files(("RF/"),full.names = TRUE)
fractions <- c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95)
Balanced <- matrices[grep("Matrix",matrices)]
Hundreds  <- matrices[grep("_100",matrices)]
x <- as.matrix(read.table(Hundreds[1]))

averages <- vector()
deviations <- vector()
averageList <- list()
for(i in 1:nrow(x)){
  good <- which(x[i,]>0)
  print(length(good))
  averages[i] <- median(as.numeric(x[i,good]))
  averageList[[i]] <- as.numeric(x[i,good])
  deviations[i] <- sd(as.numeric(x[i,good]))
}

# It says averages, but it's set to median right now.
averages1 <- averages
pdf(file='FBA_Curve.pdf')
par(mar=c(5,5,5,5))
plot(fractions,averages1,type='l',ylim=c(.5,1),
     xlab="Training Set Fraction",ylab="Balanced Accuracy",lwd=4,col='black',
     cex.main=1.5,cex.lab=2.5,cex.axis=2)
minor.tick(nx=4,ny=2,tick.ratio=.5)

saturations <- vector()
sizes <- c("_10","_20","_30","_40","_50","_60","_70","_80","_90","_100")

maxAverages <- vector()
maxSds <- vector()
curve20 <- vector()
for(j in 1:9){
  temp <-  matrices[grep(sizes[j],matrices)]
  x <- as.matrix(read.table(temp[1]))
  averages <- vector()
  deviations <- vector()
  averageList <- list()
  for(i in 1:nrow(x)){
    good <- which(x[i,]>0)
    print(length(good))
    # change median to mean to see
    # the results are pretty much the same
    averages[i] <- median(as.numeric(x[i,good])) #mean(as.numeric(x[i,good]))
    averageList[[i]] <- as.numeric(x[i,good])
    deviations[i] <- sd(as.numeric(x[i,good]))
  }
  maxAverages[j] <- median(averageList[[11]])
  maxSds[j] <- sd(averageList[[11]])
  saturation <-1- averages/averages[11]
  which(abs(saturation) < .01)
  saturations[j] <- (averages[4]-averages[1])/4
  points(fractions,averages,type='l',col=col_set[j*10],lwd=3)
  #points(deviations,type='l')
  #points(fractions[1],averages[1],,pch=paste(1:9)[j])
  #points(fractions[11],averages[11],,pch=paste(1:9)[j])
  if(j==2){
    print(averages)
    curve20 <- averages
  } 
  if(j == 1){
    curve10 <- averages
  }
}
dev.off()
pdf(file='FBA_Curve_Legend.pdf')
par(mar=c(5,20,5,8))
legend_image <- as.raster(matrix(rev(col_set), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = FALSE,xlab = '', ylab = '',cex.main=1.5)
text(x=1.5, y = seq(0,90,l=750), labels = seq(10,90,by=10),cex=1.5)
rasterImage(legend_image, 0, -.10, 1,1)
dev.off()


### Overlay the 20 member learning curve on the 14 member E coli curve
pdf(file='Ecoli_learningCurve.pdf')
par(mar=c(5,5,5,5))
plot(fractions_E,medianRFs,ylim=c(.45,1),xlab='Training Set Fraction',type='l',
     ylab='Balanced Accuracy',cex.lab=2.5,cex.main=1.5,lwd=3,cex.axis=2)
points(fractions,curve20,type='l',col=col_set[2*10],lwd=3,lty=3)
legend('bottomright',c('20 Organisms'),lty=3,lwd=3,col=col_set[20])
dev.off()
