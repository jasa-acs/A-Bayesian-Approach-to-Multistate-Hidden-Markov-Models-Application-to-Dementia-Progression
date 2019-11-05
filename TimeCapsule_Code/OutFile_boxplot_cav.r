#library(expm)
library(matrixStats)

args <- commandArgs(TRUE)
Dir <- args[1]
NumSims <- as.integer(args[2])
population <- args[3]

steps  = 15000
burnin = 10000


cmap <- matrix(c( 1, 4, 7,10,13,  16,17,18,19,20,21,
 				  2, 5, 8,11,14,   0, 0, 0, 0, 0, 0,
 				  3, 6, 9,12,15,   0, 0, 0, 0, 0, 0), 3, 11, byrow=TRUE)

parEst_hmm <- NULL
parEst_mle <- NULL
parEst_msm <- NULL
center <- NULL
meanSampleSize <- NULL
upper_hmm <- NULL
lower_hmm <- NULL
upper_mle <- NULL
lower_mle <- NULL
upper_msm <- NULL
lower_msm <- NULL
for(k in 1:NumSims){
	 
	if( file.exists( paste0(Dir,'Output_hmm',toString(k),'.rda') ) ){
		
		
		
		load(paste0(Dir,'Output_hmm',toString(k),'.rda'))
		parChains <- Output_hmm$fit$parChains
		print(Output_hmm$fit$accept)
		temp <- NULL
		temp_upper_hmm <- NULL
		temp_lower_hmm <- NULL
		for(l in 1:6){  
			param <- unlist(parChains[[l]])[burnin:steps,]
			temp <- c( temp, colMedians(param) )  
			
			for(r in 1:ncol(param)){
				temp_upper_hmm <- c( temp_upper_hmm, quantile( param[,r], probs=.975) )
				temp_lower_hmm <- c( temp_lower_hmm, quantile( param[,r], probs=.025) )
			}
		}
		parEst_hmm <- rbind( parEst_hmm, temp)	
		upper_hmm <- rbind( upper_hmm, temp_upper_hmm)
		lower_hmm <- rbind( lower_hmm, temp_lower_hmm)
	
	
		load(paste0(Dir,'Output_mle',toString(k),'.rda'))
		temp_upper_mle <- NULL
		temp_lower_mle <- NULL
		#for(l in 1:6){  			
		#	for(r in 1:ncol(param)){
		#		temp_upper_mle <- c( temp_upper_mle, quantile( param[,r], probs=.975) )
		#		temp_lower_mle <- c( temp_lower_mle, quantile( param[,r], probs=.025) )
		#	}
		#}
		parEst_mle <- rbind( parEst_mle, Output_mle$fit$coef)	
		#upper_mle <- rbind( upper_mle, temp_upper_mle)
		#lower_mle <- rbind( lower_mle, temp_lower_mle)
	
	
		load(paste0(Dir,'Output_msm',toString(k),'.rda'))
		temp <- NULL
		temp_upper_msm <- NULL
		temp_lower_msm <- NULL
		for(r in c(5,13,10,14,15)){  
			for(l in 1:3){  
				temp <- c( temp, Output_msm$Qmatrices[[l]][r])  
				temp_upper_msm <- c( temp_upper_msm, 
									 Output_msm$Qmatrices[[l]][r] + 1.96*Output_msm$QmatricesSE[[l]][r])
				temp_lower_msm <- c( temp_lower_msm, 
									 Output_msm$Qmatrices[[l]][r] - 1.96*Output_msm$QmatricesSE[[l]][r])
			}  
		}
		for(r in c(5,2,10,7)){
				temp <- c( temp, Output_msm$Ematrices[[1]][r])  
				temp_upper_msm <- c( temp_upper_msm, Output_msm$EmatricesU[[1]][r])
				temp_lower_msm <- c( temp_lower_msm, Output_msm$EmatricesL[[1]][r])
		}
		temp <- c( temp, Output_msm$opt$par[20:21])
		temp_upper_msm <- c( temp_upper_msm, Output_msm$ci[36:37,2])
		temp_lower_msm <- c( temp_lower_msm, Output_msm$ci[36:37,1])

		parEst_msm <- rbind( parEst_msm, temp)
		upper_msm <- rbind( upper_msm, temp_upper_msm)
		lower_msm <- rbind( lower_msm, temp_lower_msm)
	
		
		
		load(paste0(Dir,'meanYears',toString(k),'.rda'))
		center <- c( center, meanYears)
		load(paste0(Dir,'sampleSize',toString(k),'.rda'))
		print(N)
		meanSampleSize <- c( meanSampleSize, N)

	}
}
NumSims <- length(center)
cat('Number of output files found = ',NumSims,'\n')
cat('Mean Sample Size = ',mean(meanSampleSize),'\n')






#----------------------------------------------------------------------------------------------------------------
# Transition rate parameter boxplots ----------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


labels <- c('State 1 (well)   --->   State 2 (mild)',
			'State 1 (well)   --->   State 4 (dead)',
			'State 2 (mild)   --->   State 3 (severe)',
			'State 2 (mild)   --->   State 4 (dead)',
			'State 3 (severe)   --->   State 4 (dead)')

# betaMat is used to generate the data in Simulate_cav.r.
trueBetaMat <- matrix(c(-2.54,  0.11, -0.56,
						-2.94, -0.24,  0.15,
						-1.10, -0.15, -0.03,
						-3.92,  0.23,  0.21,
						-2.12,  0.08,  1.17), nrow=5, byrow=TRUE)
					

index <- cmap[1,1:5] # Indices from cmap.
truth <- t(trueBetaMat)[index]

if(population=='TRUE'){ pdf(paste0(Dir,'qmat_intercept_pop.pdf')) } else{ pdf(paste0(Dir,'qmat_intercept.pdf')) }
par(mfrow=c(2, 2), mar=c(5, 2.5, 2.5, 2.5))  # mar = c(bottom, left, top, right)
for(l in 1:5){
		
	intercepts_hmm <- parEst_hmm[,index[l]] - center*trueBetaMat[l,2]
	intercepts_mle <- parEst_mle[,index[l]] - center*trueBetaMat[l,2]
	intercepts_msm <- parEst_msm[,index[l]] - center*trueBetaMat[l,2]
	
	covg_hmm <- round( sum( lower_hmm[,index[l]] - center*trueBetaMat[l,2] <= truth[l] & 
			 		 truth[l] <= upper_hmm[,index[l]] - center*trueBetaMat[l,2] ) / NumSims, 3)
 	covg_msm <- round( sum( lower_msm[,index[l]] - center*trueBetaMat[l,2] <= truth[l] & 
 			 		 truth[l] <= upper_msm[,index[l]] - center*trueBetaMat[l,2] ) / NumSims, 3)
    
	if(as.logical(population)){ 
		boxplot(intercepts_hmm, intercepts_mle, intercepts_msm, names=c('Bayes','MLE','MLE (msm)'),ylab=NA, 
				main=labels[l])
	} else{
		boxplot(intercepts_hmm, intercepts_msm, names=c('Bayes','MLE (msm)'),ylab=NA, main=labels[l])
	}
	abline(h=truth[l],col='green',lwd=4,lty=5)
	mtext(paste0('.95 coverage = ',toString(covg_hmm),' Bayes, ',toString(covg_msm),' MLE (msm)'), 
		  side=1, line=2.5, cex=.75)
}
dev.off()


index <- cmap[2,1:5] # Indices from cmap.
truth <- t(trueBetaMat)[index]

if(population=='TRUE'){ pdf(paste0(Dir,'qmat_iyears_pop.pdf')) } else{ pdf(paste0(Dir,'qmat_iyears.pdf')) }
par(mfrow=c(2, 2), mar=c(5, 2.5, 2.5, 2.5))  # mar = c(bottom, left, top, right)
for(l in 1:5){
	
	covg_hmm <- round( sum( lower_hmm[,index[l]] <= truth[l] & truth[l] <= upper_hmm[,index[l]] ) / NumSims, 3)
 	covg_msm <- round( sum( lower_msm[,index[l]] <= truth[l] & truth[l] <= upper_msm[,index[l]] ) / NumSims, 3)
	
	if(as.logical(population)){ 
		boxplot(parEst_hmm[,index[l]],parEst_mle[,index[l]],parEst_msm[,index[l]],
				names=c('Bayes','MLE','MLE (msm)'), ylab=NA, main=labels[l])
	} else{
		boxplot(parEst_hmm[,index[l]],parEst_msm[,index[l]], names=c('Bayes','MLE (msm)'),ylab=NA,main=labels[l])
	}
	abline(h=truth[l],col='green',lwd=4,lty=5)
	mtext(paste0('.95 coverage = ',toString(covg_hmm),' Bayes, ',toString(covg_msm),' MLE (msm)'), 
		  side=1, line=2.5, cex=.75)
}
dev.off()


index <- cmap[3,1:5] # Indices from cmap.
truth <- t(trueBetaMat)[index]

if(population=='TRUE'){ pdf(paste0(Dir,'qmat_sex_pop.pdf')) } else{ pdf(paste0(Dir,'qmat_sex.pdf')) }
par(mfrow=c(2, 2), mar=c(5, 2.5, 2.5, 2.5))  # mar = c(bottom, left, top, right)
for(l in 1:5){
	
	covg_hmm <- round( sum( lower_hmm[,index[l]] <= truth[l] & truth[l] <= upper_hmm[,index[l]] ) / NumSims, 3)
 	covg_msm <- round( sum( lower_msm[,index[l]] <= truth[l] & truth[l] <= upper_msm[,index[l]] ) / NumSims, 3)
	
	if(as.logical(population)){ 
		boxplot(parEst_hmm[,index[l]],parEst_mle[,index[l]],parEst_msm[,index[l]],
				names=c('Bayes','MLE','MLE (msm)'), ylab=NA, main=labels[l])
	} else{
		boxplot(parEst_hmm[,index[l]],parEst_msm[,index[l]], names=c('Bayes','MLE (msm)'),ylab=NA,main=labels[l])
	}
	abline(h=truth[l],col='green',lwd=4,lty=5)
	mtext(paste0('.95 coverage = ',toString(covg_hmm),' Bayes, ',toString(covg_msm),' MLE (msm)'), 
		  side=1, line=2.5, cex=.75)
}
dev.off()


#----------------------------------------------------------------------------------------------------------------
# Misclassification parameter boxplots --------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


trueErrors <- c( .01, .24, .06, .11)

labels <- c('P( observed state 2 | true state 1 )',
			'P( observed state 1 | true state 2 )',
			'P( observed state 3 | true state 2 )',
			'P( observed state 2 | true state 3 )')

index <- cmap[1,6:9] # Indices from cmap.

# logistic transforms of the parameter etimates.
vals_hmm <- matrix(0,NumSims,length(index))
vals_mle <- matrix(0,NumSims,length(index))
vals_msm <- matrix(0,NumSims,length(index))
for(k in 1:NumSims){
	vals_hmm[k,] <- c( exp(parEst_hmm[k,index[1]]) / sum( c( 1, exp(parEst_hmm[k,index[1]])) ),
					   exp(parEst_hmm[k,index[2:3]]) / sum( c( 1, exp(parEst_hmm[k,index[2:3]])) ),
					   exp(parEst_hmm[k,index[4]]) / sum( c( 1, exp(parEst_hmm[k,index[4]])) ) )
				
   	vals_mle[k,] <- c( exp(parEst_mle[k,index[1]]) / sum( c( 1, exp(parEst_mle[k,index[1]])) ),
   					   exp(parEst_mle[k,index[2:3]]) / sum( c( 1, exp(parEst_mle[k,index[2:3]])) ),
   					   exp(parEst_mle[k,index[4]]) / sum( c( 1, exp(parEst_mle[k,index[4]])) ) )
			
	vals_msm[k,] <- c( exp(parEst_msm[k,index[1]]) / sum( c( 1, exp(parEst_msm[k,index[1]])) ),
					   exp(parEst_msm[k,index[2:3]]) / sum( c( 1, exp(parEst_msm[k,index[2:3]])) ),
					   exp(parEst_msm[k,index[4]]) / sum( c( 1, exp(parEst_msm[k,index[4]])) ) )
}

truth <- c( log(trueErrors[1] / .99),
		    log(trueErrors[2:3] / .70),
		    log(trueErrors[4] / .89) )

if(population=='TRUE'){ pdf(paste0(Dir,'misclass_pop.pdf')) } else{ pdf(paste0(Dir,'misclass.pdf')) }
par(mfrow=c(2, 2), mar=c(5, 2.5, 2.5, 2.5))  # mar = c(bottom, left, top, right)
for(l in 1:4){
	
	
	covg_hmm <- round( sum( lower_hmm[,index[l]] <= truth[l] & truth[l] <= upper_hmm[,index[l]] ) / NumSims, 3)
 	covg_msm <- round( sum( lower_msm[,index[l]] <= trueErrors[l] & 
							trueErrors[l] <= upper_msm[,index[l]] ) / NumSims, 3)
	
	if(as.logical(population)){
		boxplot(vals_hmm[,l],vals_mle[,l],vals_msm[,l],names=c('Bayes','MLE','MLE (msm)'),ylab=NA,main=labels[l])
	} else{
		boxplot(vals_hmm[,l],vals_msm[,l],names=c('Bayes','MLE (msm)'), ylab=NA, main=labels[l])
	}
	abline(h=trueErrors[l],col='green',lwd=4,lty=5)
	mtext(paste0('.95 coverage = ',toString(covg_hmm),' Bayes, ',toString(covg_msm),' MLE (msm)'), 
		  side=1, line=2.5, cex=.75)
}
dev.off()


#----------------------------------------------------------------------------------------------------------------
# Initial probability parameter boxplots ------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


trueInitProbs <- c( .04, .01)

labels <- c('P( initial state 2 )','P( initial state 3 )')

index <- cmap[1,10:11] # Indices from cmap.

# logistic transforms of the parameter etimates.
vals_hmm <- matrix(0,NumSims,length(index))
vals_mle <- matrix(0,NumSims,length(index))
vals_msm <- matrix(0,NumSims,length(index))
for(k in 1:NumSims){
	vals_hmm[k,] <- exp(parEst_hmm[k,index]) / sum( c( 1, exp(parEst_hmm[k,index]), 0) )
	vals_mle[k,] <- exp(parEst_mle[k,index]) / sum( c( 1, exp(parEst_mle[k,index]), 0) )
	vals_msm[k,] <- exp(parEst_msm[k,index]) / sum( c( 1, exp(parEst_msm[k,index]), 0) )
}

truth <- log(trueInitProbs / .95)

if(population=='TRUE'){ pdf(paste0(Dir,'initialP_pop.pdf')) } else{ pdf(paste0(Dir,'initialP.pdf')) }
par(mfrow=c(2, 2), mar=c(5, 2.5, 2.5, 2.5))  # mar = c(bottom, left, top, right)
for(l in 1:2){
	
	
	covg_hmm <- round( sum( lower_hmm[,index[l]] <= truth[l] & truth[l] <= upper_hmm[,index[l]] ) / NumSims, 3)
 	covg_msm <- round( sum( lower_msm[,index[l]] <= trueInitProbs[l] & 
					 trueInitProbs[l] <= upper_msm[,index[l]] ) / NumSims, 3)
	
	if(as.logical(population)){
		boxplot(vals_hmm[,l],vals_mle[,l],vals_msm[,l],names=c('Bayes','MLE','MLE (msm)'),ylab=NA,main=labels[l])
	} else{
		boxplot(vals_hmm[,l],vals_msm[,l],names=c('Bayes','MLE (msm)'), ylab=NA, main=labels[l])
	}
	abline(h=trueInitProbs[l],col='green',lwd=4,lty=5)
	mtext(paste0('.95 coverage = ',toString(covg_hmm),' Bayes, ',toString(covg_msm),' MLE (msm)'), 
		  side=1, line=2.5, cex=.75)
}
dev.off()



