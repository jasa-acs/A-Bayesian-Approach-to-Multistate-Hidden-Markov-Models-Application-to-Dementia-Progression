library(expm)
library(splines)

args <- commandArgs(TRUE)
Dir <- args[1]
NumSeeds <- as.integer(args[2])
scale <- args[3]
seeds <- 1:NumSeeds

steps  <- 30000
burnin <- 25000

states <- c("A-N-", "A+N-", "A-N+", "A+N+", "A-Dem", "A+Dem", "Dead")
plotCol <- c('palegreen3','plum','skyblue','purple1','purple4','firebrick3')

cmap <- matrix(c( 1,12,17,22,27,32,37,42,42,42,42,58,63,  68,74,75,76,77,78,79,     80:87,  88,89,90,
 				  2,13,18,23,28,33,38,43,43,43,43,59,64,  69,69,69,69,69,69, 0,  rep(0,8),   0, 0, 0,
 				  3,14,19,24,29,34,39,44,44,44,44,60,65,  70,70,70,70,70,70, 0,  rep(0,8),   0, 0, 0,
 				  4,15,20,25,30,35,40,45,45,45,45,61,66,  71,71,71,71,71,71, 0,  rep(0,8),   0, 0, 0,
 				  5,16,21,26,31,36,41,46,46,46,46,62,67,  72,72,72,72,72,72, 0,  rep(0,8),   0, 0, 0,
				  6, 0, 0, 0, 0, 0, 0,47,47,47,47, 0, 0,  73,73,73,73,73,73, 0,  rep(0,8),   0, 0, 0,
 				  7, 0, 0, 0, 0, 0, 0,48,48,48,48, 0, 0,   0, 0, 0, 0, 0, 0, 0,  rep(0,8),   0, 0, 0,
 				  8, 0, 0, 0, 0, 0, 0,49,49,49,49, 0, 0,   0, 0, 0, 0, 0, 0, 0,  rep(0,8),   0, 0, 0,
 				  9, 0, 0, 0, 0, 0, 0,50,50,50,50, 0, 0,   0, 0, 0, 0, 0, 0, 0,  rep(0,8),   0, 0, 0,
 				 10, 0, 0, 0, 0, 0, 0,51,51,51,51, 0, 0,   0, 0, 0, 0, 0, 0, 0,  rep(0,8),   0, 0, 0,
 				 11, 0, 0, 0, 0, 0, 0,52,52,52,52, 0, 0,   0, 0, 0, 0, 0, 0, 0,  rep(0,8),   0, 0, 0,
 				  0, 0, 0, 0, 0, 0, 0,53,53,53,53, 0, 0,   0, 0, 0, 0, 0, 0, 0,  rep(0,8),   0, 0, 0,
 				  0, 0, 0, 0, 0, 0, 0,54,54,54,54, 0, 0,   0, 0, 0, 0, 0, 0, 0,  rep(0,8),   0, 0, 0,
 				  0, 0, 0, 0, 0, 0, 0,55,55,55,55, 0, 0,   0, 0, 0, 0, 0, 0, 0,  rep(0,8),   0, 0, 0,
 				  0, 0, 0, 0, 0, 0, 0,56,56,56,56, 0, 0,   0, 0, 0, 0, 0, 0, 0,  rep(0,8),   0, 0, 0,
				  0, 0, 0, 0, 0, 0, 0,57,57,57,57, 0, 0,   0, 0, 0, 0, 0, 0, 0,  rep(0,8),   0, 0, 0
				  ), 16, 31, byrow=TRUE)
				 
parGroups_Index = list()
parGroups_Index[[1]] = cmap[,c(1,4,7)][cmap[,c(1,4,7)]>0]
parGroups_Index[[2]] = cmap[,c(2,3)][cmap[,c(2,3)]>0]
parGroups_Index[[3]] = cmap[,c(5,6)][cmap[,c(5,6)]>0]
parGroups_Index[[4]] = cmap[,8][cmap[,8]>0][1:7]
parGroups_Index[[5]] = cmap[,12:13][cmap[,12:13]>0]
parGroups_Index[[6]] = unique(cmap[,14:20][cmap[,14:20]>0])
parGroups_Index[[7]] = cmap[,21:31][cmap[,21:31]>0]

chainList <- NULL
for(k in 1:NumSeeds){
	
	load(paste0(Dir,'Output',toString(seeds[k]),'.rda'))
	
	chain <- matrix(0, steps-burnin+1, 90)
	for(l in c(1,2,3, 5,6,7)){  chain[,parGroups_Index[[l]]] <- Output$fit$parChains[[l]][burnin:steps,]  }
	# Add zeros for non-dem to death terms.
	chain[,cmap[,8][cmap[,8]>0]] <- cbind( Output$fit$parChains[[4]][burnin:steps,], matrix(0,steps-burnin+1,9))
		
	chainList[[k]] <- chain	
}
StackedChains <- do.call(rbind,chainList)





Basis_amyl <- bs( 50:120, knots=c(55,65,75,90), degree=3, intercept=TRUE)
colnames(Basis_amyl) <- c('c1_amyl','c2_amyl','c3_amyl','c4_amyl','c5_amyl','c6_amyl','c7_amyl','c8_amyl')
				 

Q <- function(iage,male,educ,apoe4,par){
	
	rateParIndex = matrix(c(NA,12,17,22,27,32,37,42,58,63,  
		 				    NA,13,18,23,28,33,38,43,59,64,  
		 				    NA,14,19,24,29,34,39,44,60,65,  
		 				    NA,15,20,25,30,35,40,45,61,66,  
		 				    NA,16,21,26,31,36,41,46,62,67), 5, 10, byrow=TRUE)
	
	q1  = exp( c(male,educ,apoe4,Basis_amyl[iage+68-49,]) %*% par[1:11] )  # Transition from state 1 to state 2.
	q2  = exp( c(1,iage,male,educ,apoe4) %*% par[rateParIndex[,2]] )  # Transition from state 1 to state 3.
	q3  = exp( c(1,iage,male,educ,apoe4) %*% par[rateParIndex[,3]] )  # Transition from state 2 to state 4.
	q4  = exp( c(1,iage,male,educ,apoe4) %*% par[rateParIndex[,4]] )  # Transition from state 3 to state 4.
	q5  = exp( c(1,iage,male,educ,apoe4) %*% par[rateParIndex[,5]] )  # Transition from state 3 to state 5.
	q6  = exp( c(1,iage,male,educ,apoe4) %*% par[rateParIndex[,6]] )  # Transition from state 4 to state 6.
	q7  = exp( c(1,iage,male,educ,apoe4) %*% par[rateParIndex[,7]] )  # Transition from state 5 to state 6.
	q8  = exp( c(1,iage,male,educ,apoe4) %*% par[rateParIndex[,8]] )  # Transition from state 1 to state 7.
	
	# Don't add the death bias here; the bias only applies to someone in the study.
	
	q9  = q8      									  				  # Transition from state 2 to state 7.
	q10 = q8   										 				  # Transition from state 3 to state 7.
	q11 = q8 										 				  # Transition from state 4 to state 7.
	q12 = exp( c(1,iage,male,educ,apoe4) %*% par[rateParIndex[,9]] )  # Transition from state 5 to state 7.
	q13 = exp( c(1,iage,male,educ,apoe4) %*% par[rateParIndex[,10]] ) # Transition from state 6 to state 7.

	qmat = matrix(c(  0, q1, q2,  0,  0,  0, q8,
					  0,  0,  0, q3,  0,  0, q9,
					  0,  0,  0, q4, q5,  0,q10,
				  	  0,  0,  0,  0,  0, q6,q11,
				  	  0,  0,  0,  0,  0, q7,q12,
				  	  0,  0,  0,  0,  0,  0,q13,
					  0,  0,  0,  0,  0,  0,  0), ncol=7, byrow=TRUE) 
	diag(qmat) = -rowSums(qmat)
	
	states = c("A-N-", "A+N-", "A-N+", "A+N+", "A-Dem", "A+Dem", "Dead")
	dimnames(qmat) = list(states, states)
	
	return(qmat)
}





#----------------------------------------------------------------------------------------------------------------
# Log death rate.
#----------------------------------------------------------------------------------------------------------------

library(survival)
DRates = survexp.mn[,,'2000']
ageDR = as.numeric(rownames(DRates))
maleDR = DRates[ageDR >= 50,'male']*365
femaleDR = DRates[ageDR >= 50,'female']*365
ageDR = ageDR[ageDR >= 50]

print(summary(lm( log(maleDR) ~ 1 + ageDR )))
print(summary(lm( log(femaleDR) ~ 1 + ageDR )))

pdf(paste0(Dir,'logDeathRates_',scale,'.pdf'))
plot( ageDR, log(femaleDR), type='l', lwd=3, col='black', xlab='Age', main='Log death rate',
	  cex.axis=1.25, cex.main=2, cex.lab=1.5, ylab='Annualized death rate (not on log-scale)', yaxt='n')
lines(ageDR,log(maleDR),type='l',lty=2, lwd=3,col='black')
axis(2, at=log(femaleDR), labels=round( exp(log(femaleDR)), 4), cex.axis=1)
dev.off()

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
# Death bias
#----------------------------------------------------------------------------------------------------------------

# Plot the b-spline using the sample means of the control points.
C <- -colMeans(StackedChains)[47]
D <- colMeans(StackedChains)[48] + C
aVec <- seq(-C,by=D,length=11)
aVec[aVec>0] <- 0
B <- aVec

# Empirically compute credible regions for the spline, sampling from the MCMC sample points.
bspline_sampleSize <- 1000
B_sample <- matrix(0,bspline_sampleSize,11)
for(k in 1:bspline_sampleSize){
	
	sample_index <- sample( 1:dim(StackedChains)[1], size=1)
	
	C <- -StackedChains[sample_index,47]
	D <- StackedChains[sample_index,48] + C
	aVec <- seq(-C,by=D,length=11)
	aVec[aVec>0] <- 0
	B_sample[k,] <- aVec
}

upper_bspline <- rep(0,11)
lower_bspline <- rep(0,11)
for(k in 1:11){
	upper_bspline[k] <- quantile( B_sample[,k], probs=.975)
	lower_bspline[k] <- quantile( B_sample[,k], probs=.025)
}

x <- 0:10
pdf(paste0(Dir,'deathBiasSpline_',scale,'.pdf'))
plot( x, lower_bspline, xlab='Years enrolled', ylab='Proportion of population death rate', main='Death bias', 
	  type='l', col=plotCol[6], panel.first=grid( NULL, NULL, col="gray", lty="dotted", lwd=3, equilogs=TRUE), 
	  yaxt='n', ylim=c(log(.25),0), cex.axis=1.25, cex.main=2, cex.lab=1.5, lwd=3, lty=2)
lines( x, upper_bspline, col=plotCol[6], lwd=3, lty=2)
lines( x, B, col=plotCol[6], lwd=3, lty=1)
axis(2, at=c(lower_bspline[1],B), labels=round(exp(c(lower_bspline[1],B)),2), cex.axis=1)
dev.off()

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
# Cubic spline for transition from A-N- (state 1) to A+N- (state 2)
#----------------------------------------------------------------------------------------------------------------

# Plot the b-spline using the sample means of the control points.
cubicSpline <- Basis_amyl %*% colMeans(StackedChains)[4:11]

# Empirically compute credible regions for the spline, sampling from the MCMC sample points.
bspline_sampleSize = 1000
cubicSpline_sample = matrix(0,bspline_sampleSize,71)
for(k in 1:bspline_sampleSize){
	
	sample_index <- sample( 1:dim(StackedChains)[1], size=1)

	cubicSpline_sample[k,] = Basis_amyl %*% StackedChains[sample_index,4:11]
}

upper_cubicSpline = rep(0,71)
lower_cubicSpline = rep(0,71)
for(k in 1:71){
	upper_cubicSpline[k] <- quantile( cubicSpline_sample[,k], probs=.975)
	lower_cubicSpline[k] <- quantile( cubicSpline_sample[,k], probs=.025)
}

x <- 50:100
pdf(paste0(Dir,'cubicSpline_s1_to_s2_',scale,'.pdf'))
plot( x, lower_cubicSpline[1:51], xlab='Age', ylab='log rate', main='State 1 (A-N-)   --->   State 2 (A+N-)', 
	  type='l', col=plotCol[6], panel.first=grid( NULL, NULL, col="gray", lty="dotted", lwd=3, equilogs=TRUE), 
	  ylim=c(-8.5,-2.2), cex.axis=1.25, cex.main=2, cex.lab=1.5, lwd=3, lty=2)
lines( x, upper_cubicSpline[1:51], col=plotCol[6], lwd=3, lty=2)
lines( x, cubicSpline[1:51], col=plotCol[6], lwd=3, lty=1)
dev.off()

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
# Resonse data histograms
#----------------------------------------------------------------------------------------------------------------
library(mclust)
load('synthetic_MCSA_data.rda')

parEst <- colMeans(StackedChains)


pdf(paste0(Dir,'density_MMSE_',scale,'.pdf'))
mmse <- useData$mmse[!is.na(useData$mmse)]

hold <- Mclust( mmse, G=6, modelNames='E')$parameters

# Must un-normalize the mmse data.
density1 <- hold$pro[6]*dnorm(seq(15,32,by=.01),mean=2.170923*parEst[68]+27.68659,sd=2.170923*exp(parEst[79])^.5)
density2 <- hold$pro[5]*dnorm(seq(15,32,by=.01),mean=2.170923*parEst[74]+27.68659,sd=2.170923*exp(parEst[79])^.5)
density3 <- hold$pro[4]*dnorm(seq(15,32,by=.01),mean=2.170923*parEst[75]+27.68659,sd=2.170923*exp(parEst[79])^.5)
density4 <- hold$pro[3]*dnorm(seq(15,32,by=.01),mean=2.170923*parEst[76]+27.68659,sd=2.170923*exp(parEst[79])^.5)
density5 <- hold$pro[2]*dnorm(seq(15,32,by=.01),mean=2.170923*parEst[77]+27.68659,sd=2.170923*exp(parEst[79])^.5)
density6 <- hold$pro[1]*dnorm(seq(15,32,by=.01),mean=2.170923*parEst[78]+27.68659,sd=2.170923*exp(parEst[79])^.5)

hist(2.170923*mmse+27.68659, freq=FALSE, xlim=c(15,32),main='Distribution of observed MMSE', xlab=NA, ylab=NA,
	 col='seashell2', border='gray', breaks=30, cex.axis=1.25, cex.main=2)

lines( seq(15,32,by=.01), density1, col=plotCol[6], lwd=3)
lines( seq(15,32,by=.01), density2, col=plotCol[5], lwd=3)
lines( seq(15,32,by=.01), density3, col=plotCol[4], lwd=3)
lines( seq(15,32,by=.01), density4, col=plotCol[3], lwd=3)
lines( seq(15,32,by=.01), density5, col=plotCol[2], lwd=3)
lines( seq(15,32,by=.01), density6, col=plotCol[1], lwd=3)

lines( seq(15,32,by=.01), density1 + density2 + density3 + density4 + density5 +density6, col=4, lwd=2, lty=2)
legend('topleft', legend=c('state 1','state 2','state 3','state 4','state 5','state 6','normal mixture density'),
	    col=c(plotCol[6:1],4), lwd=3, lty=c(1,1,1,1,1,1,2), cex=1.25, box.lty=0)
dev.off()



pdf(paste0(Dir,'density_pib_',scale,'.pdf'))
pib <- exp(useData$lpib[!is.na(useData$lpib)])+1

hold <- Mclust( log(pib-1), G=2, modelNames='E')$parameters

density1 <- hold$pro[1]*dlnorm( seq(1,3.5,by=.01)-1, mean=parEst[80], sd=exp(parEst[82])^.5)
density2 <- hold$pro[2]*dlnorm( seq(1,3.5,by=.01)-1, mean=parEst[81], sd=exp(parEst[82])^.5)

hist(pib, freq=FALSE, xlim=c(1,3.5), main='Observed PIB density estimates', xlab=NA, ylab=NA,
	 col='seashell2', border='gray', breaks=sqrt(length(pib)), cex.axis=1.25, cex.main=2, xaxt='n')

lines( seq(1,3.5,by=.01), density1, col=plotCol[1], lwd=3)
mean <- exp(parEst[80])+1
upper_sd <- exp(parEst[80] + exp(parEst[82])^.5)+1
lower_sd <- exp(parEst[80] - exp(parEst[82])^.5)+1
segments(x0=mean, y0=0, x1=mean, 
		 y1=hold$pro[1]*dlnorm( mean-1, mean=parEst[80], sd=exp(parEst[82])^.5),col=plotCol[1],lwd=6,lty=5)
segments(x0=upper_sd, y0=0, x1=upper_sd, 
	     y1=hold$pro[1]*dlnorm( upper_sd-1, mean=parEst[80], sd=exp(parEst[82])^.5),col=plotCol[1],lwd=3,lty=5)
segments(x0=lower_sd, y0=0, x1=lower_sd, 
	     y1=hold$pro[1]*dlnorm( lower_sd-1, mean=parEst[80], sd=exp(parEst[82])^.5),col=plotCol[1],lwd=3,lty=5)
xaxisVals <- c(lower_sd,mean,upper_sd)

lines( seq(1,3.5,by=.01), density2, col=plotCol[2], lwd=3)
mean <- exp(parEst[81])+1
upper_sd <- exp(parEst[81] + exp(parEst[82])^.5)+1
lower_sd <- exp(parEst[81] - exp(parEst[82])^.5)+1
segments(x0=mean, y0=0, x1=mean, 
	 	 y1=hold$pro[2]*dlnorm( mean-1, mean=parEst[81], sd=exp(parEst[82])^.5),col=plotCol[2],lwd=6,lty=5)
segments(x0=upper_sd, y0=0, x1=upper_sd, 
     	 y1=hold$pro[2]*dlnorm( upper_sd-1, mean=parEst[81], sd=exp(parEst[82])^.5),col=plotCol[2],lwd=3,lty=5)
segments(x0=lower_sd, y0=0, x1=lower_sd, 
     	 y1=hold$pro[2]*dlnorm( lower_sd-1, mean=parEst[81], sd=exp(parEst[82])^.5),col=plotCol[2],lwd=3,lty=5)
xaxisVals <- c(xaxisVals, lower_sd,mean,upper_sd)

axis(1, at=xaxisVals, labels=round(xaxisVals,2), cex.axis=1.25, las=2)

lines( seq(1,3.5,by=.01), density1 + density2, col=4, lwd=2, lty=2)
legend('topright', legend=c('low burden','high burden','normal mixture density'),
	    col=c(plotCol[1:2],4), lwd=3, lty=c(1,1,2), cex=1.25, box.lty=0)
dev.off()



pdf(paste0(Dir,'density_thickness_',scale,'.pdf'))
thickness <- useData$thickness[!is.na(useData$thickness)]

hold <- Mclust( thickness, G=2, modelNames='E')$parameters

density1 <- hold$pro[1]*dnorm( seq(1.75,3.4,by=.01), mean=parEst[84], sd=exp(parEst[85])^.5)
density2 <- hold$pro[2]*dnorm( seq(1.75,3.4,by=.01), mean=parEst[83], sd=exp(parEst[85])^.5)

hist(thickness, freq=FALSE, xlim=c(1.75,3.4), main='Observed Thickness density estimates', xlab=NA, ylab=NA,
	 col='seashell2', border='gray', breaks=sqrt(length(thickness)), cex.axis=1.25, cex.main=2, xaxt='n')

lines( seq(1.75,3.4,by=.01), density2, col=plotCol[1], lwd=3)
mean <- parEst[83]
upper_sd <- parEst[83] + exp(parEst[85])^.5
lower_sd <- parEst[83] - exp(parEst[85])^.5
segments(x0=mean, y0=0, x1=mean, 
	 	 y1=hold$pro[2]*dnorm( mean, mean=parEst[83], sd=exp(parEst[85])^.5),col=plotCol[1],lwd=6,lty=5)
segments(x0=upper_sd, y0=0, x1=upper_sd, 
     	 y1=hold$pro[2]*dnorm( upper_sd, mean=parEst[83], sd=exp(parEst[85])^.5),col=plotCol[1],lwd=3,lty=5)
segments(x0=lower_sd, y0=0, x1=lower_sd, 
     	 y1=hold$pro[2]*dnorm( lower_sd, mean=parEst[83], sd=exp(parEst[85])^.5),col=plotCol[1],lwd=3,lty=5)
xaxisVals <- c(lower_sd,mean,upper_sd)


lines( seq(1.75,3.4,by=.01), density1, col=plotCol[2], lwd=3)
mean <- parEst[84]
upper_sd <- parEst[84] + exp(parEst[85])^.5
lower_sd <- parEst[84] - exp(parEst[85])^.5
segments(x0=mean, y0=0, x1=mean, 
		 y1=hold$pro[1]*dnorm( mean, mean=parEst[84], sd=exp(parEst[85])^.5),col=plotCol[2],lwd=6,lty=5)
segments(x0=upper_sd, y0=0, x1=upper_sd, 
	     y1=hold$pro[1]*dnorm( upper_sd, mean=parEst[84], sd=exp(parEst[85])^.5),col=plotCol[2],lwd=3,lty=5)
segments(x0=lower_sd, y0=0, x1=lower_sd, 
	     y1=hold$pro[1]*dnorm( lower_sd, mean=parEst[84], sd=exp(parEst[85])^.5),col=plotCol[2],lwd=3,lty=5)
xaxisVals <- c(xaxisVals, lower_sd,mean,upper_sd)

axis(1, at=unique(round(xaxisVals,2)), labels=unique(round(xaxisVals,2)), cex.axis=1.25, las=2)

lines( seq(1.75,3.4,by=.01), density1 + density2, col=4, lwd=2, lty=2)
legend('topleft', legend=c('low burden','high burden','normal mixture density'),
	    col=c(plotCol[1:2],4), lwd=3, lty=c(1,1,2), cex=1.25, box.lty=0)
dev.off()

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
# Dementia diagnosis misclassification estimates
#----------------------------------------------------------------------------------------------------------------

eta <- colMeans(StackedChains)[86:87]

DemMisclass <- matrix(0,2,2)
rownames(DemMisclass) <- c('NonDem','Dem')
colnames(DemMisclass) <- c('diagnosed NonDem','diagnosed Dem')
DemMisclass[1,] <- c( 1, exp(eta[1])) / (1 + exp(eta[1]))
DemMisclass[2,] <- c( exp(eta[2]), 1) / (1 + exp(eta[2]))
print(DemMisclass)

CredInt_upper <- matrix(0,2,2)
rownames(CredInt_upper) <- c('NonDem','Dem')
colnames(CredInt_upper) <- c('diagnosed NonDem','diagnosed Dem')
CredInt_lower <- matrix(0,2,2)
rownames(CredInt_lower) <- c('NonDem','Dem')
colnames(CredInt_lower) <- c('diagnosed NonDem','diagnosed Dem')

upper <- quantile( StackedChains[,86], probs=.975)
lower <- quantile( StackedChains[,86], probs=.025)
CredInt_upper[1,] <- c( 1, exp(upper)) / (1 + exp(upper))
CredInt_lower[1,] <- c( 1, exp(lower)) / (1 + exp(lower))

upper <- quantile( StackedChains[,87], probs=.975)
lower <- quantile( StackedChains[,87], probs=.025)
CredInt_upper[2,] <- c( exp(upper), 1) / (1 + exp(upper))
CredInt_lower[2,] <- c( exp(lower), 1) / (1 + exp(lower))

print(CredInt_upper)
print(CredInt_lower)
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
# Probability evolution plots
#----------------------------------------------------------------------------------------------------------------

ProbEvol_dementia <- function(state,male,educ,apoe4,par){
	
	probEvo = rep(0,50)
	
	# Compute the probability of transitioning to state 5 given state 1 at age=50 and not dead.
	P = c(1,0,0,0,0,0,0)
	for(iage in 50:99) {
		P = P %*% expm( Q(iage-68,male,0,apoe4,par) )
		probEvo[iage-49] = sum(P[state]) / sum(P[1:6])
	}
	
	return(probEvo)
}


ProbEvol_alzheimer <- function(male,educ,apoe4,par){
	
	dt = 1/12
	probEvo = rep(0,50)
	for(t in 51:100){  # For book keeping, it is assumed here that t is an integer. 
		
		cummulProb = 0
		# h is the exact tranistion time to state 4.  Must leave at least one instance of time, dt, to
		# transition to state 6 by time t (hence, h cannot exceed t-dt). 
		for(h in seq( from=50+dt, to=t-dt, by=dt)){
			
			# Compute the probability of transitioning from state 1 to states 1-3 from age 50 to h.		
			part = unique( c( 50:floor(h), seq( floor(h), h, by=dt)) )
			dt_part = diff(part)
			part = setdiff(part,floor(h))
			
			P = c(1,0,0,0,0,0,0)
			for( k in 1:length(part)) P = P %*% expm( dt_part[k] * Q( floor(part[k])-68, male,educ,apoe4,par) )
			# Finally, compute the probability of transitioning to state 4 exactly at time h.
			prob_s4 = P[1:3] %*% Q( floor(h)-68, male,educ,apoe4,par)[1:3,4] * dt
			
			
			# Compute the probability of transitioning to state 6 given state 4 at age=h and not dead.
			part = unique( c( seq( h, ceiling(h), by=dt), ceiling(h):t) )
			dt_part = diff(part)
			part = setdiff(part,t)
			
			P = c(0,0,0,1,0,0,0)
			for( k in 1:length(part)) P = P %*% expm( dt_part[k] * Q( floor(part[k])-68, male,educ,apoe4,par) )
			prob_s6_given_s4 = P[6] / sum(P[1:6])
			
			# Sum over discretization of all possible exact times of trans to state 4, denoted by the h values.
			cummulProb = cummulProb + prob_s6_given_s4 * prob_s4
		}
		probEvo[t-50] = cummulProb
		print(t)
	}	
	
	return(probEvo)
}



vals <- rbind( c(0,0,0,colMeans(StackedChains)),
			   c(1,0,0,colMeans(StackedChains)),
			   c(0,0,1,colMeans(StackedChains)),
			   c(1,0,1,colMeans(StackedChains)))
				   
labels <- c('female, apeo4 negative',
			'male, apeo4 negative',
			'female, apeo4 positive',
			'male, apeo4 positive')

pdf(paste0(Dir,'probEvol_',scale,'.pdf'))
par(mfrow=c(2, 2), mar=c(5, 2.5, 2.5, 2.5))  # mar = c(bottom, left, top, right)
for(k in 1:4){
	plot( NULL, ylab='Probability (given not dead)', xlab='Age', ylim=c(0,.3), xlim=c(50,100), 
		  main=labels[k], panel.first=grid( NULL, NULL, col="gray", lty="dotted", lwd=3, equilogs=TRUE))
	
	Dementia <- ProbEvol_dementia( state=5:6, male=vals[k,1], educ=vals[k,2], apoe4=vals[k,3], par=vals[k,4:90])
	AminusDem <- ProbEvol_dementia( state=5, male=vals[k,1], educ=vals[k,2], apoe4=vals[k,3], par=vals[k,4:90])
	alzheimer <- ProbEvol_alzheimer( male=vals[k,1], educ=vals[k,2], apoe4=vals[k,3], par=vals[k,4:90])

	lines( 51:100, Dementia, col=plotCol[3], lwd=3)
	lines( 51:100, AminusDem, col=plotCol[4], lwd=3)
	lines( 51:100, alzheimer, col=plotCol[6], lwd=3)
	if(k==1) legend( 'topleft', c('Dementia (any form)','A-Dem (state 5)',"Alzheimer's"), 
					 col=plotCol[c(3,4,6)], bty='n', lwd=3)
}
dev.off()
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
# State space heat maps
#----------------------------------------------------------------------------------------------------------------
library(diagram)
library(animation)



position <- matrix( c(.13,.61,
				      .37,.85,
					  .37,.37,
					  .61,.61,
					  .61,.13,
					  .85,.37,
					  .85,.85), byrow=TRUE, ncol=2)

vals <- rbind( c(0,0,0,colMeans(StackedChains)),
			   c(1,0,0,colMeans(StackedChains)),
			   c(0,0,1,colMeans(StackedChains)),
			   c(1,0,1,colMeans(StackedChains)))
				   
labels <- c(' year old female, apeo4 negative',
			' year old male, apeo4 negative',
			' year old female, apeo4 positive',
			' year old male, apeo4 positive')
			

maxHeatCol = 350
topoColors <- colorRampPalette(c('red','yellow','green','blue'), bias=2.5)
topoColors <- paste0(topoColors(maxHeatCol), "80")

pdf(paste0(Dir,'heatMaps_',scale,'.pdf'))
par(mfrow=c(2, 2), mar=c(0,0,0,0) + 0.8, oma = c(0,0,0,0) + 0.8)  # mar = c(bottom, left, top, right)
for(t in 50:90){	  
	for(k in 1:4){
	  				  
		qmat <- Q( iage=floor(t)-68, male=vals[k,1], educ=vals[k,2], apoe4=vals[k,3], par=vals[k,4:90])
		qmat = round( qmat, 3)*(qmat>0)
		qmat[c(1,3,4),7] <- 0
		qmat[qmat>0] <- round( 1/qmat[qmat>0], 1)

		qmatCol <- matrix(0,7,7)
		for(i in 1:7)  for(j in 1:7)  if(round(qmat[i,j],0) > 0)  qmatCol[i,j] <- topoColors[round(qmat[i,j],0)]

		plotmat( t(qmat), pos=position, curve=0, name=states, lwd=1, box.lwd=2, cex.txt=1, box.type="ellipse", 
				 box.prop=0.5, arr.type="triangle", arr.pos=.75, arr.lcol=t(qmatCol), arr.col=t(qmatCol), 
				 arr.lwd=5, arr.width=.5, shadow.size=0, endhead=TRUE, main=paste0( toString(t), labels[k]))
	}
	print(t)
}
dev.off()



pdf(paste0(Dir,'heatLegend_',scale,'.pdf'))
par(pin = c(6.5,.25))
plot( seq( 1, maxHeatCol, by=1), rep( 0, maxHeatCol), col=topoColors, pch=19, cex=1.5, 
	  xlab=NA, ylab=NA, yaxt='n')
dev.off()
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------




