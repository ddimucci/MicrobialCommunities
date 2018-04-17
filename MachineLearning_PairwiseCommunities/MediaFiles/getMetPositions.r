# THIS FILE GOES THROUGH THE FLUX LOGS OF EVERY INTERACTION AND IDENTIFIES WHICH METABOLITES WERE
# COMPETED FOR, COOPERATED FOR ETC.

# Feb 16 2018
# When we rank the metabolites, how far will we have to go down the list before we encounter the first two mechanisms?

# MAKE SURE TO MODIFY EACH OF THESE PATHS TO POINT TO WHERE YOU STORED THE NAMED FILE
Elements <- read.table("/projectnb/cometsfba/dimucci/Pairwise/Features/Curated_Element_Names.txt")
Elements <- Elements[,1]
Elements <- as.vector(Elements)
JointElements <- c(Elements,Elements)
JointElements[2084:4166] <- paste0(JointElements[2084:4166],'.1')
use <- grep('EX',JointElements)
Exchanges <- JointElements[use]
Interactions <- read.table('/projectnb/cometsfba/dimucci/Pairwise/Data/RYTable.txt')
Data <- read.table('/projectnb/cometsfba/dimucci/Pairwise/Features/ConcatenatedVectors.txt')

Monos <- which(Interactions[,1]==Interactions[,2])
exData <- Data[-Monos,use]
colnames(exData) <- Exchanges
Interactions <- Interactions[-Monos,]

# I CREATED TWO FEATURE CONTRIBUTION MATRICES
# FCs_full is the feature contribution matrix obtained from the full 9900 element dataset
# FCs_tiny is the matrix obtained from the balanced data set built to classify the 420 positive samples 
X <- read.table('/projectnb/cometsfba/dimucci/Pairwise/Features/FCs_full.txt')
Y <- read.table('/projectnb/cometsfba/dimucci/Pairwise/Features/FCs_tiny.txt')
colnames(Y) <- colnames(X)

Positions<- matrix(0,nrow=9900,ncol=14)
colnames(Positions) <- c('FirstNeg','SecondNeg','FirstPos','SecondPos','TotalMechs','TotalNeg','TotalPos','FalseMech1','FalseMech2',
        'FalseMech3','FalseMech4',"FluxRank1",'FluxRank2','FluxRankCombine')

CompetitionIndex <- rep(0,194)
CrossFeedIndex <- rep(0,194)
CompeteMatrix <- matrix(0,ncol=194,nrow=9900)
FeedMatrix <- CompeteMatrix


MetNames <- matrix(0,nrow=9900,ncol=4)
colnames(MetNames) <- c("FirstNeg","SecondNeg","FirstPos","SecondPos")

# Keep track of all positions. Might as well again.
Neg_positions <- matrix(0,nrow=9900,ncol=194)
Pos_positions <- Neg_positions
Potential_positions <- Neg_positions

Neg_Indices <- Neg_positions
Pos_Indices <- Neg_positions

Neg_Flux <- Neg_positions
Pos_Flux <- Neg_positions

# Move to the directory where the output is kept
setwd("/projectnb/cometsfba/dimucci/Pairwise/Output")

Directories <- list.files(full.name=T)

# For each directory move into it and read the flux logs
# Use the flux logs to identify which metabolites are the DIRECT mechanisms of competition or cooperation
count = 0

for(j in 1:length(Directories)){
        files <- list.files(Directories[j])
        # temp.txt indicates that GSEA has been run on this interaction. Which means that the mechanisms are found
        if('temp.txt' %in% files){
                setwd(Directories[j])
                models <- read.csv('/projectnb/cometsfba/dimucci/Pairwise/Models/modelCodes.txt',header=F)
                fluxes <- list.files()
                pair <- fluxes[grep('\\.mat',fluxes)]
                pair <- gsub('_xml.mat','',pair)
                pair[1] <- grep(pair[1],models[,1])
                pair[2] <- grep(pair[2],models[,1])

                # If there are 2 models in here proceed
                if(length(pair) > 1){
                        fluxes <- fluxes[grep('species',fluxes)]

                        if(length(fluxes > 1)){
                                flux1 <- read.csv('species1.csv')
                                flux2 <- read.csv('species2.csv')
                        }


                        colnames(flux1)=''
                        colnames(flux2)=''

                        fluxes <- read.csv('parsedFluxes.txt')

                        target1 <- which(Interactions[,1]==pair[1] & Interactions[,2]==pair[2])
                        target2 <- which(Interactions[,1]==pair[2] & Interactions[,2]==pair[1])


                        # Map fluxes to compound IDs
                        both <- rbind(flux1,flux2)
                        C = paste0(both[1:nrow(flux1),1],'')
                        L=paste0(both[(nrow(flux1)+1):nrow(both),1],'.1')
                        both <- c(C,L)
                        fluxes <- cbind(both,fluxes)

                        # Identify exchanges both models have
                        keep = which(both%in%Exchanges)

                        # Fluxes of the immediately relevant predictors
                        predictors <- fluxes[keep,]

                        second <- grep('\\.1',predictors[,1])
                        second <- which(predictors[,5]==2)
                        secondFlux <- predictors[second,]
                        firstFlux <-predictors[-second,]


                        # Find the reactions that both organisms have, these are the most likely sources 
                        # of the interaction
                        firstP <- which(exData[target1,1:194]==1)
                        secondP <- which(exData[target1,195:388]==1)
                        bothPresent <- intersect(firstP,secondP) # potential mechanisms

                        firstP1 <- which(exData[target2,1:194]==1)
                        secondP1 <- which(exData[target2,195:388]==1)
                        bothPresent1 <- intersect(firstP,secondP) # potential mechanisms

                        # Find the reactions that both had and keep them
                        # These are the potential competition mechanisms
                        common <- rep(0,nrow(firstFlux))
                        keep <- rep(0,nrow(firstFlux))
                        for(i in 1:nrow(firstFlux)){
                          if(length(grep(firstFlux[i,1],secondFlux[,1]))>0){
                            common[i]=grep(firstFlux[i,1],secondFlux[,1])
                            keep[i]=1
                          }
                        }

                        # Arrange them side by side, these are the metabolites both organisms have transporters for.
                        # Interactions have to happen through at least one of these pairs
                        combineFluxes <- cbind(firstFlux[which(keep==1),],secondFlux[common,])
                        netFluxes <- combineFluxes[,7]+combineFluxes[,14]

                        # Define a bunch of indices here. 
                        # We really only care about competeIndex and feedIndex
                        # But record the rest of the interaction types just in case
                        # find the reactions where cross feeding occurred
                        crossFeed <- which(combineFluxes[,7]*combineFluxes[,14]<0)
                        feedIndex <-which(Exchanges %in% combineFluxes[crossFeed,1])
                        crossFeedNames <- combineFluxes[crossFeed,1]

                        fluxNames <- combineFluxes[,1]
                        # find the reactions where competition occurred 
                        competitors <- which(combineFluxes[,7]<0 & combineFluxes[,14]<0)
                        competeIndex <- which(Exchanges %in% combineFluxes[competitors,1])
                        competeNames <- combineFluxes[competitors,1]
                        names(competitors) <- competeNames
                        nonCompete <- bothPresent[-which(bothPresent %in% competeIndex)]

                        # Record the competitive metabolites
                        if(length(competeIndex) > 0){
                                CompeteMatrix[target1,competeIndex] = 1
                                CompeteMatrix[target2,competeIndex] = 1
                        }
                        if(length(feedIndex) > 0){
                                FeedMatrix[target1,feedIndex] = 1
                                FeedMatrix[target2,feedIndex] = 1
                        }

                        # First eater
                        firstEat <- which(combineFluxes[,7] < 0 & combineFluxes[,14]==0)
                        secondEat <- which(combineFluxes[,7] == 0 & combineFluxes[,14] <0)
                        firstSecrete <- which(combineFluxes[,7] > 0 & combineFluxes[,14]==0)
                        secondSecrete <- which(combineFluxes[,7] == 0 & combineFluxes[,14]>0)
                        neitherCare <- which(combineFluxes[,7] == 0 & combineFluxes[,14] ==0)

                        # Record how many of each thing there was
                        Positions[target1,5] <- length(bothPresent)
                        Positions[target1,6] <- length(competeIndex)
                        Positions[target1,7] <- length(feedIndex)

                        Positions[target2,5] <- length(bothPresent)
                        Positions[target2,6] <- length(competeIndex)
                        Positions[target2,7] <- length(feedIndex)

                        # Record the indices of the competitive interactions and cross-feeding so we 
                        # can see which things are commonly fought over 
                        if(length(competeIndex) > 0){
                                Neg_Indices[target1,1:length(competeIndex)] <- competeIndex
                                Neg_Indices[target2,1:length(competeIndex)] <- competeIndex
                                competeOrder <- vector()
                                for(p in 1:length(competeIndex)){
                                        competeOrder[p] <- which(combineFluxes[competitors,1] == Exchanges[competeIndex[p]])
                                }

                                Neg_Flux[target1,competeIndex] <- combineFluxes[competitors,7][competeOrder]
                                Neg_Flux[target2,competeIndex] <- combineFluxes[competitors,14][competeOrder]
                        }

                        # Do the same for cross feed
                        if(length(feedIndex) > 0){
                                Pos_Indices[target1,1:length(feedIndex)] <- feedIndex
                                Pos_Indices[target2,1:length(feedIndex)] <- feedIndex
                                feedOrder <- vector()
                                }
                                Pos_Flux[target1,feedIndex] <- combineFluxes[crossFeed,7][feedOrder]
                                Pos_Flux[target2,feedIndex] <- combineFluxes[crossFeed,14][feedOrder]
                        }
                        # Now that we know the indices of both mechanism types we can rank the features.



                        #Find the positions we would encounter the competitor mechanisms at if we were to rank the metabolites
                        # by predictive influence and go down it.
                        colnames(X) <- Exchanges
                        colnames(Y) <- Exchanges
                        competeNames <- colnames(X)[competeIndex]
                        feedNames <- colnames(X)[feedIndex]
                        bothNames <- colnames(X)[bothPresent]

                        # Find the rank of competitive metabolites
                        if(length(competeIndex)> 0){
                                if(Interactions[target1,3] <= 0){
                                        FCs <- X[target1,1:194]+X[target1,195:388]
                                        Ranked <- sort(FCs)
                                        CompetitionIndex[competeIndex] <- CompetitionIndex[competeIndex] + 1
                                } else {
                                        target11 = which(rownames(Y)==rownames(Interactions)[target1])
                                        FCs <- Y[target11,1:194]+Y[target11,195:388]
                                        Ranked <- sort(FCs)
                                }
                                bothPos <- which(colnames(Ranked) %in% bothNames)
                                Potential_positions[target1,1:length(bothPos)] <- bothPos
                                Potential_positions[target2,1:length(bothPos)] <- bothPos

                                positions <- which(colnames(Ranked) %in% competeNames)
                                Neg_positions[target1,1:length(positions)]=positions

                                #print(positions)
                                #print(c('TARGET!',target1))
                                Positions[target1,1] <- positions[1]
                                Positions[target1,2] <- positions[2]

                                MetNames[target1,1] <- colnames(Ranked)[positions[1]]
                                MetNames[target1,2] <- colnames(Ranked)[positions[2]]
        
                                ##### RECORD THE RANKS OF THE TOP PREDICTOR IN TARGET1
                                fluxMet <- which(combineFluxes[competitors,1] == colnames(Ranked)[positions[1]])
                                competeRank1 <- rank(combineFluxes[competitors,7])[fluxMet]
                                competeRank2 <- rank(combineFluxes[competitors,14])[fluxMet]
                                competeRankCombine <- rank(combineFluxes[competitors,7]+combineFluxes[competitors,14])[fluxMet]
                                
                                #CompetitionIndex[competeIndex] <- CompetitionIndex[competeIndex] + 1
                                  
                                Positions[target1,8] <- length(which(bothPos < Positions[target1,1]))
                                Positions[target1,9] <- length(which(bothPos < Positions[target1,2]))
                                Positions[target1,12] <- competeRank1
                                Positions[target1,13] <- competeRank2
                                Positions[target1,14] <- competeRankCombine
                                
                                if(Interactions[target2,3] <= 0){
                                        FCs <- X[target2,1:194]+X[target2,195:388]
                                        Ranked <- sort(FCs)
                                        CompetitionIndex[competeIndex] <- CompetitionIndex[competeIndex] + 1
                                } else {
                                        target22 = which(rownames(Y)==rownames(Interactions)[target2])
                                        FCs <- Y[target22,1:194]+Y[target22,195:388]
                                        Ranked <- sort(FCs)
                                }
                                Positions[target2,5] <- length(bothPresent)
                                bothPos <- which(colnames(Ranked) %in% bothNames)
                                positions <- which(colnames(Ranked) %in% competeNames)
                                Neg_positions[target2,1:length(positions)]=positions

                                Positions[target2,1] <- positions[1]
                                Positions[target2,2] <- positions[2]

                                MetNames[target2,1] <- colnames(Ranked)[positions[1]]
                                MetNames[target2,2] <- colnames(Ranked)[positions[2]]

                                Positions[target2,8] <- length(which(bothPos < Positions[target2,1]))
                                Positions[target2,9] <- length(which(bothPos < Positions[target2,2]))
                                
                                ##### RECORD THE RANKS OF THE TOP PREDICTOR IN TARGET2
                                fluxMet <- which(combineFluxes[competitors,1] == colnames(Ranked)[positions[1]])
                                competeRank1 <- rank(combineFluxes[competitors,7])[fluxMet]
                                competeRank2 <- rank(combineFluxes[competitors,14])[fluxMet]
                                competeRankCombine <- rank(combineFluxes[competitors,7]+combineFluxes[competitors,14])[fluxMet]
                                
                                Positions[target2,8] <- length(which(bothPos < Positions[target2,1]))
                                Positions[target2,9] <- length(which(bothPos < Positions[target2,2]))
                                Positions[target2,12] <- competeRank2 
                                Positions[target2,13] <- competeRank1
                                Positions[target2,14] <- competeRankCombine
                        }

                        # Find the rank of cross-fed metabolites
                        if(length(feedIndex)> 0){
                                if(Interactions[target1,3] <= 0){
                                        FCs <- X[target1,1:194]+X[target1,195:388]
                                        Ranked <- sort(FCs)
                                } else {
                                        count <- count + 1
                                        #print(count)
                                        target11 = which(rownames(Y)==rownames(Interactions)[target1])
                                        FCs <- Y[target11,1:194]+Y[target11,195:388]
                                        Ranked <- sort(FCs)
                                        #print(colnames(rev(Ranked))[1])
                                }
                                bothPos <- which(colnames(Ranked) %in% bothNames)
                                Potential_positions[target1,1:length(bothPos)] <- bothPos
                                Potential_positions[target2,1:length(bothPos)] <- bothPos
                                positions <- which(colnames(Ranked) %in% feedNames)
                                Pos_positions[target1,1:length(positions)]=positions

                                Positions[target1,3] <- max(positions)
                                MetNames[target1,3] <- colnames(rev(Ranked))[1]#colnames(Ranked)[positions[which.max(positions)]]
                                Positions[target1,10] <- length(which(bothPos > Positions[target1,3]))
                                
                                ##### RECORD THE RANKS OF THE TOP PREDICTOR IN TARGET1
                                fluxMet <- which(combineFluxes[crossFeed,1] == colnames(Ranked)[positions[1]])
                                crossFeedRank1 <- rank(combineFluxes[crossFeed,7])[fluxMet]
                                crossFeedRank2 <- rank(combineFluxes[crossFeed,14])[fluxMet]
                                crossFeedRankCombine <- rank(combineFluxes[crossFeed,7]+combineFluxes[crossFeed,14])[fluxMet]
                                
                                CrossFeedIndex[feedIndex] <- CrossFeedIndex[feedIndex] + 1
                                
                                Positions[target2,8] <- length(which(bothPos < Positions[target1,1]))
                                Positions[target2,9] <- length(which(bothPos < Positions[target1,2]))
                                Positions[target2,12] <- crossFeedRank1 
                                Positions[target2,13] <- crossFeedRank2
                                Positions[target2,14] <- crossFeedRankCombine
                                if(length(positions) > 1){
                                        Positions[target1,4] <- max(positions[-which.max(positions)])
                                        
                                        MetNames[target1,4] <- colnames(rev(Ranked))[2]#colnames(Ranked)[positions[which.max(positions[-which.max(positions)])]]
                                        Positions[target1,11] <- length(which(bothPos > Positions[target1,4]))
                                }
                                

                                if(Interactions[target2,3] <= 0){
                                        FCs <- X[target2,1:194]+X[target2,195:388]
                                        Ranked <- sort(FCs)
                                } else {
                                        target22 = which(rownames(Y)==rownames(Interactions)[target2])
                                        FCs <- Y[target22,1:194]+Y[target22,195:388]
                                        Ranked <- sort(FCs)
                                }
                                bothPos <- which(colnames(Ranked) %in% bothNames)
                                positions <- which(colnames(Ranked) %in% feedNames)
                                Pos_positions[target2,1:length(positions)]=positions
                                Positions[target2,3] <- max(positions)
                                MetNames[target2,3] <- colnames(rev(Ranked))[1]
                                Positions[target2,10] <- length(which(bothPos > Positions[target2,3]))
                                ##### RECORD THE RANKS OF THE TOP PREDICTOR IN TARGET2
                                fluxMet <- which(combineFluxes[crossFeed,1] == colnames(Ranked)[positions[1]])
                                crossFeedRank1 <- rank(combineFluxes[crossFeed,7])[fluxMet]
                                crossFeedRank2 <- rank(combineFluxes[crossFeed,14])[fluxMet]
                                crossFeedRankCombine <- rank(combineFluxes[crossFeed,7]+combineFluxes[crossFeed,14])[fluxMet]
                                
                                CrossFeedIndex[feedIndex] <- CrossFeedIndex[feedIndex] + 1
                                
                                Positions[target2,8] <- length(which(bothPos < Positions[target2,1]))
                                Positions[target2,9] <- length(which(bothPos < Positions[target2,2]))
                                Positions[target2,12] <- crossFeedRank2 # Relative to its own fluxes
                                Positions[target2,13] <- crossFeedRank1 
                                Positions[target2,14] <- crossFeedRankCombine
                                if(length(positions) > 1){
                                        Positions[target2,4] <- max(positions[-which.max(positions)])

                                        MetNames[target2,4] <- colnames(rev(Ranked))[2]
                                        Positions[target2,11] <- length(which(bothPos > Positions[target2,4]))
                                }
                        }
                }
        }
        setwd("/projectnb/cometsfba/dimucci/Pairwise/Output")
        print(j)
}

UsageIndex <- rbind(CompetitionIndex,CrossFeedIndex)
colnames(UsageIndex) <- colnames(FCs)
# MAKE SURE TO RENAME THESE PATHS TOO
write.csv(UsageIndex, file='/projectnb/cometsfba/dimucci/Pairwise/Results/RANK_RESULTS/UsageIndex.csv')
write.csv(Positions,file='/projectnb/cometsfba/dimucci/Pairwise/Results/RANK_RESULTS/Positions_mechs.csv')
write.csv(CompeteMatrix,file='/projectnb/cometsfba/dimucci/Pairwise/Results/RANK_RESULTS/CompeteMatrix.csv')
write.csv(FeedMatrix,file='/projectnb/cometsfba/dimucci/Pairwise/Results/RANK_RESULTS/FeedMatrix.csv')
write.csv(Neg_Flux,file='/projectnb/cometsfba/dimucci/Pairwise/Results/RANK_RESULTS/NegFlux.csv')
write.csv(Pos_Flux,file='/projectnb/cometsfba/dimucci/Pairwise/Results/RANK_RESULTS/PosFlux.csv')
