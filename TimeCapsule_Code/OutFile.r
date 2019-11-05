#library(splines)


args <- commandArgs(TRUE)
Dir <- args[1]
NumSeeds <- as.integer(args[2])
scale <- args[3]

steps  <- 30000
burnin <- 25000

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
	chain[,cmap[,8][cmap[,8]>0]] <- cbind( Output$fit$parChains[[4]][burnin:steps,], matrix(0,steps-burnin+1,9))
		
	chainList[[k]] <- chain	
}
StackedChains <- do.call(rbind,chainList)
print(colMeans(StackedChains[,4:11]))
# Compute death bias slope.
StackedChains[,48] <- StackedChains[,48] - StackedChains[,47]


filesNames <- c(paste0(Dir,'trace_qmat_intercept_',scale,'.pdf'),
				paste0(Dir,'trace_qmat_iage_',scale,'.pdf'),
				paste0(Dir,'trace_qmat_male_',scale,'.pdf'),
				paste0(Dir,'trace_qmat_educ_',scale,'.pdf'),
				paste0(Dir,'trace_qmat_apoe4_',scale,'.pdf'),
				paste0(Dir,'trace_mmse_',scale,'.pdf'),
				paste0(Dir,'trace_Resp_',scale,'.pdf'),
				paste0(Dir,'trace_initialProbs_',scale,'.pdf'))

indexList <- list(c(12,17,22,27,32,37,42,58,63,47,48),
				  c(4:11,13,18,23,28,33,38,43,59,64),
				  c( 1,14,19,24,29,34,39,44,60,65),
				  c( 2,15,20,25,30,35,40,45,61,66),
				  c( 3,16,21,26,31,36,41,46,62,67),
			      c(68,74,75,76,77,78,69:73,79),
				  c(80:87),
				  c(88:90))

transLabels <- c('State 1 (A-N-)   --->   State 2 (A+N-)',
			     'State 1 (A-N-)   --->   State 3 (A-N+)',
			     'State 2 (A+N-)   --->   State 4 (A+N+)',
			     'State 3 (A-N+)   --->   State 4 (A+N+)',
			     'State 3 (A-N+)   --->   State 5 (DemA-)',
		         'State 4 (A+N+)   --->   State 6 (DemA+)',
			     'State 5 (DemA-)   --->   State 6 (DemA+)',
			     'States 1-4 (NonDem)   --->   State 7 (Dead)',
			     'State 5 (DemA-)   --->   State 7 (Dead)',
			     'State 6 (DemA+)   --->   State 7 (Dead)')
				 
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
				 				  
				

				
for(r in 1:8){
	
	pdf(filesNames[r])
	index <- indexList[[r]]
	labels <- labelsList[[r]]
	par(mfrow=c(3, 2))
	for(l in 1:length(index)){
		
		plot( NULL, ylab=NA, xlab=NA, xlim=c(1,steps-burnin+1), ylim=range(StackedChains[,index[l]]),
			 main=labels[l])
		
		for(k in 1:NumSeeds){
			if(r == 1 & l == 11){
				lines( chainList[[k]][,index[11]] - chainList[[k]][,index[10]], type='l', col=k)		
			} else{
				lines( chainList[[k]][,index[l]], type='l', col=k)  
			}
		}
		
		parMean = round( mean(StackedChains[,index[l]]), 4)
		parMedian = round( median(StackedChains[,index[l]]), 4)
		upper = quantile( StackedChains[,index[l]], prob=.975)
		lower = quantile( StackedChains[,index[l]], prob=.025)
		
		hist( StackedChains[,index[l]], breaks=sqrt(as.integer(args[2])*(steps-burnin)), ylab=NA, main=NA,
			  freq=FALSE, xlab=paste0('Mean =',toString(parMean),' Median =',toString(parMedian)))
		abline( v=upper, col='red', lwd=2, lty=2)
		abline( v=lower, col='purple', lwd=2, lty=2)
	}
	dev.off()

}





