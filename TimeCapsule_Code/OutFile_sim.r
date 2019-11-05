#library(splines)

args <- commandArgs(TRUE)
Dir <- args[1]
NumSeeds <- as.integer(args[2])

steps  <- 15000
burnin <- 10000
seeds <- 1:NumSeeds

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


trueBetaMat <- matrix(c(          NA,		   NA,-0.331664059, 0.014273922, 0.945136344,
						-5.000000000, 0.094236594, 0.141479641,-0.095807126, 0.486230191,
						-4.200000000, 0.093838248,-0.116771101,-0.098580487, 0.399355964,
						-3.881605509, 0.065857892,-0.330382263, 0.019152541, 0.444423684,
						-4.148148694, 0.075589995, 0.248030366, 0.033380640, 0.765503215,
						-2.473865801, 0.045208611, 0.082975736, 0.054272656, 0.256422760,
						-2.739931535, 0.047242201, 0.050003054, 0.033386975, 0.239889316,
						-4.329458534, 0.110780340, 0.387452143,-0.042883267,-0.008324821,
						-3.500000000, 0.120000000, 0.360000000,-0.007003006, 0.104718992,
						-1.100000000, 0.060000000, 0.224925174, 0.003157103, 0.236539477), ncol=5,byrow=TRUE)


coef <- c(-0.014554096,-0.089207557, 0.058688151,-0.004627584, 0.057826362)
intercept <- c( 0.538414931, 0.178435764,-0.672421219,-1.827065032,-2.341332296,-5.208797471)
truthMMSE <- c( intercept, coef, -1.261838927)

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
	print(Output$fit$accept)
	
	chain <- matrix(0, steps-burnin+1, 90)
	for(l in c(1,2,3, 5,6,7)){  chain[,parGroups_Index[[l]]] <- Output$fit$parChains[[l]][burnin:steps,]  }
	# Add zeros for non-dem to death terms.
	chain[,cmap[,8][cmap[,8]>0]] <- cbind(Output$fit$parChains[[4]][burnin:steps,1:6],
										  Output$fit$parChains[[4]][burnin:steps,7]-
										  Output$fit$parChains[[4]][burnin:steps,6],matrix(0,steps-burnin+1,9))

	chainList[[k]] <- chain	
}
StackedChains <- do.call(rbind,chainList)



filesNames <- c(paste0(Dir,'trace_qmat_intercept.pdf'),
				paste0(Dir,'trace_qmat_iage.pdf'),
				paste0(Dir,'trace_qmat_male.pdf'),
				paste0(Dir,'trace_qmat_educ.pdf'),
				paste0(Dir,'trace_qmat_apoe4.pdf'),
				paste0(Dir,'trace_mmse.pdf'),
				paste0(Dir,'trace_Resp.pdf'),
				paste0(Dir,'trace_initialProbs.pdf'))
				
filesNames_boxplot <- c(paste0(Dir,'boxplot_qmat_intercept.pdf'),
						paste0(Dir,'boxplot_qmat_iage.pdf'),
						paste0(Dir,'boxplot_qmat_male.pdf'),
						paste0(Dir,'boxplot_qmat_educ.pdf'),
						paste0(Dir,'boxplot_qmat_apoe4.pdf'),
						paste0(Dir,'boxplot_mmse.pdf'),
						paste0(Dir,'boxplot_Resp.pdf'),
						paste0(Dir,'boxplot_initialProbs.pdf'))

indexList <- list(c(12,17,22,27,32,37,42,58,63,47,48),
				  c(4:11,13,18,23,28,33,38,43,59,64),
				  c( 1,14,19,24,29,34,39,44,60,65),
				  c( 2,15,20,25,30,35,40,45,61,66),
				  c( 3,16,21,26,31,36,41,46,62,67),
			      c(68,74,75,76,77,78,69:73,79),
				  c(80:87),
				  c(88:90))

transLabels <- c('(A-N-) ---> (A+N-)',
			     '(A-N-) ---> (A-N+)',
			     '(A+N-) ---> (A+N+)',
			     '(A-N+) ---> (A+N+)',
			     '(A-N+) ---> (DemA-)',
		         '(A+N+) ---> (DemA+)',
			     '(DemA-) ---> (DemA+)',
			     '(NonDem) ---> (Dead)',
			     '(DemA-) ---> (Dead)',
			     '(DemA+) ---> (Dead)')
				 
labelsList <- list(	c(transLabels[2:10], 'death bias baseline','death bias slope'), 
					c('c1_amyl','c2_amyl','c3_amyl','c4_amyl','c5_amyl','c6_amyl','c7_amyl','c8_amyl', 
					  transLabels[2:10]), 
					transLabels, transLabels, transLabels, 

					c('State 1 (A-N-)',
					  'State 2 (A+N-)',
					  'State 3 (A-N+)',
					  'State 4 (A+N+)',
					  'State 5 (DemA-)',
					  'State 6 (DemA+)',
					  'age',
					  'male',
					  'educ',
					  'apoe4',
					  'ntests',
					  'log(sigma^2)'),

					c('Mean PIB | A-',
		  			  'Mean PIB | A+',
		  			  'log(sigma^2) for PIB',
		  			  'Mean Thickness | N-',
		  			  'Mean Thickness | N+',
		  			  'log(sigma^2) for Thickness',
		  			  'logit P( observed Dem | NonDem )',
		  			  'logit P( observed NonDem | Dem )'),

					c('mlogit P( initial state 2 )',
					  'mlogit P( initial state 3 )',
					  'mlogit P( initial state 4 )'))

 					 

truthList <- list( c(trueBetaMat[2:10,1], -1.242275027, -1.112106935+1.242275027),
                   c(-6.580919931,-6.501250885,-6.427228628,-4.673031945,
					 -3.965086053,-3.111756265,-1.950950409, 3.588549810, trueBetaMat[2:10,2]),
				   trueBetaMat[,3],
				   trueBetaMat[,4],
				   trueBetaMat[,5],
				   truthMMSE,
				   c(-1.084010116, 0.110684536,-2.240498159,
 					  2.655601472, 2.498948044,-3.869343130,
 					 -4.712799832,-2.163822374),
				   c(-3.9835,-4.9499,-6.5185))




# Compute the limits of 95 percent credible intervals for each parameter, in each chain, and posterior means.
upper <- matrix(0,NumSeeds,max(cmap))
lower <- matrix(0,NumSeeds,max(cmap))
postMeans <- matrix(0,NumSeeds,max(cmap))
for(k in 1:NumSeeds){
for(p in 1:max(cmap)){
	upper[k,p] <- quantile( chainList[[k]][,p], probs=.975)
	lower[k,p] <- quantile( chainList[[k]][,p], probs=.025)
}

postMeans[k,] <- colMeans(chainList[[k]])
}





#----------------------------------------------------------------------------------------------------------------
# Produce trace plots of the stacked MCMC chains, after convergence
#----------------------------------------------------------------------------------------------------------------
for(r in 1:8){
	
	pdf(filesNames[r])
	index <- indexList[[r]]
	labels <- labelsList[[r]]
	truth <- truthList[[r]]
	par(mfrow=c(3, 2))
	for(l in 1:length(index)){
		
		plotRange <- c( min( min(StackedChains[,index[l]]), truth[l]), 
						max( max(StackedChains[,index[l]]), truth[l]))
						
		plot( NULL, ylab=NA, xlab=NA, xlim=c(1,steps-burnin+1), ylim=range(StackedChains[,index[l]]),
			 main=labels[l])
			 
		for(k in 1:NumSeeds) lines( chainList[[k]][,index[l]], type='l', col=k)  
		
	
		parMean = round( mean(StackedChains[,index[l]]), 4)
		parMedian = round( median(StackedChains[,index[l]]), 4)
		
		hist( StackedChains[,index[l]], breaks=sqrt(as.integer(NumSeeds)*(steps-burnin)), ylab=NA, main=NA,
			  freq=FALSE, xlab=paste0('Mean =',toString(parMean),' Median =',toString(parMedian)),
			  xlim=plotRange)
		abline( v=truth[l], col='green', lwd=3, lty=1)
	}
	dev.off()

}
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------



				
				
#----------------------------------------------------------------------------------------------------------------
# Produce boxplots with 95 percent coverage proportions
#----------------------------------------------------------------------------------------------------------------
for(r in 1:8){	
	
	pdf(filesNames_boxplot[r])
	index <- indexList[[r]]
	labels <- labelsList[[r]]
	truth <- truthList[[r]]
	par(mfrow=c(2, 2), mar=c(5, 2.5, 2.5, 2.5))  # mar = c(bottom, left, top, right)
	for(l in 1:length(index)){
		
		covg <- sum( lower[,index[l]] <= truth[l] & truth[l] <= upper[,index[l]] ) / NumSeeds
		
		plotRange <- c( min( min(postMeans[,index[l]]), truth[l]), 
						max( max(postMeans[,index[l]]), truth[l]))
		
		boxplot( postMeans[,index[l]], names=c('Posterior means'), ylim=plotRange, ylab=NA, main=labels[l], 
				 outline=TRUE, cex.axis=1.25, cex.main=1.5)
		
		abline(h=truth[l],col='green',lwd=4,lty=5)
		mtext(paste0('.95 coverage = ',toString(round(covg,3))), side=1, line=2.5, cex=1.25)
	}
	dev.off()

}
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------



