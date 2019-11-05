library(MASS)
library(mvtnorm)

		 
mcmcRoutine <- function(par, fn, cmap, prior.means, prior.std, steps, burnin){
	
    # 'par' contains the initial coefficient values, which are also the prior means if 'prior.means' is missing.
    npar = length(par)
	
    if(missing(prior.means)) prior.means <- par
    if(length(prior.means) != npar) stop("wrong length for prior.means")

    if(missing(cmap)) stop("cmap argument is required")
    if(max(cmap) != npar) stop("cmap does not agree with the number of parameters")
		
    if(missing(prior.std)) stop("prior.std argument is required")
    if(length(prior.std) == 1) prior.std <- rep(prior.std, npar)
    if(length(prior.std) != npar) stop("wrong length for prior.std")
	
	
	# Partition the paramenters to be updated in groups -----------------------------------------------------
	# The group index.
	parGroups = 1:6
	
	# The positive components of 'cmap' correspond to the parameters to be updatded.  The positive components also correspond to the indices of 'par'.
	parGroups_Index = list()
	parGroups_Index[[1]] = cmap[,1][cmap[,1]>0]
	parGroups_Index[[2]] = cmap[,2][cmap[,2]>0]
	parGroups_Index[[3]] = cmap[,3][cmap[,3]>0]
	parGroups_Index[[4]] = cmap[,4][cmap[,4]>0]
	parGroups_Index[[5]] = cmap[,5][cmap[,5]>0]
	parGroups_Index[[6]] = cmap[,6:11][cmap[,6:11]>0]
	#------------------------------------------------------------------------------------------------------- 


	parChains = list() # This list stores the MCMC parameter vector chain for each parameter group.  
	for(l in parGroups){ parChains[[l]] = matrix(0,steps,length(parGroups_Index[[l]])) }
	
	Sigma = list()  # This list stores the proposal covariance matrices for each parameter group.
	for(l in parGroups){ Sigma[[l]] = diag( rep(1,length(parGroups_Index[[l]])) ) }

	tau = rep(.001,length(parGroups)) # These are used to tune the proposals during the burnin period.
		
	prev_logPosteriorDens = -Inf  # Initialize.
	proposal = par # Initialize.

	### The algorithm...
	accept = rep(0,length(parGroups))
	for(ttt in 1:steps){

		### Update each of the beta vectors, in turn.
		for(l in parGroups){ 
			
			# The indices of 'par' corresponding to transition 'l'.
			tempInd = parGroups_Index[[l]]
			proposal[tempInd] = par[tempInd] + mvrnorm(1,rep(0,length(tempInd)),tau[l]*Sigma[[l]])
			
			# Compute the log of the posterior density of the proposed parameter vector.
			logPosteriorDens = fn(proposal) + dmvnorm(proposal, prior.means, diag(prior.std^2), log=TRUE)
			
			# Accept/reject proposed parameters.
			logMH_Ratio = logPosteriorDens - prev_logPosteriorDens
			if( log(runif(1,0,1)) < logMH_Ratio ){ # If TRUE, accept the proposed parameter.
				
				# The updates if the proposal is accepted.
				prev_logPosteriorDens = logPosteriorDens
				par = proposal
				accept[l] = accept[l] + 1 
		
			} else{  proposal = par  }
			
			parChains[[l]][ttt,] = par[tempInd]
			
			
			# Proposal tuning scheme ----------------------------------------------------------------------------
			if(ttt < burnin){
				# During the burnin period, update the proposal covariance in each step to capture the relationships within the parameters vectors for each transition.  This helps with mixing.
				if(ttt == 100){  tau[l] = 1  }
				
				if(100 <= ttt & ttt <= 2000){  
					chain = parChains[[l]][1:ttt,]
					Sigma[[l]] = cov(chain[ !duplicated(chain),, drop=FALSE])
					
				} else if(2000 < ttt){  
					chain = parChains[[l]][(ttt-2000):ttt,]
					Sigma[[l]] = cov(chain[ !duplicated(chain),, drop=FALSE])
				}
				if( sum( is.na(Sigma[[l]]) ) > 0 ){  Sigma[[l]] = diag( length(parGroups_Index[[l]]) )  }
			
				# Tune the proposal covariance for each transition to achieve reasonable acceptance ratios.
				if(ttt %% 30 == 0){ 
					if(ttt %% 480 == 0){  
						accept[l] = 0  
		
					} else if( accept[l] / (ttt %% 480) < 0.45 ){ 
						tau[l] = (.5^2)*tau[l] 
						
					} else if( accept[l] / (ttt %% 480) > 0.55 ){ 
						tau[l] = (1.5^2)*tau[l] 
					} 
				}
			}
			#----------------------------------------------------------------------------------------------------
			

			#cat('groupNum',l,'\n')
		}
		
		# Restart the acceptance ratio at burnin.
		if(ttt == burnin){  accept = rep(0,length(parGroups)) }
		cat('------------>',ttt,'\n')
	}
	
	cat('Number of steps = ',ttt,'\n')
	accept = accept / (steps-burnin)
	print(accept)
	
	list(loglik=prev_logPosteriorDens, tau=tau, accept=accept, par=par, parGroups_Index=parGroups_Index, 
			parChains=parChains, Sigma=Sigma)
}





