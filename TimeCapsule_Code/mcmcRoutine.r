library(MASS)
library(mvtnorm)

		 
mcmcRoutine <- function(par, fn, cmap, splineParam, prior.means, prior.std, steps, burnin, designMat, 
																				   cubicSplineMat, splineMat){
	
    # 'par' contains the initial coefficient values, which are also the prior means if 'prior.means' is missing.
    npar = length(par)

    if(missing(prior.means)) prior.means <- par
    if(length(prior.means) != npar) stop("wrong length for prior.means")

    if(missing(cmap)) stop("cmap argument is required")
    if(max(cmap) != npar) stop("cmap does not agree with the number of parameters")
		
    if(missing(prior.std)) stop("prior.std argument is required")
    if(length(prior.std) == 1) prior.std <- rep(prior.std, npar)
    if(length(prior.std) != npar) stop("wrong length for prior.std")
	
	
	# Partition the paramenters to be updated in groups ---------------------------------------------------------
	# The group index.
	parGroups = 1:7
	
	# The positive components of 'cmap' correspond to the parameters to be updatded.  The positive components also correspond to the indices of 'par'.
	parGroups_Index = list()
	parGroups_Index[[1]] = cmap[,c(1,4,7)][cmap[,c(1,4,7)]>0]
	parGroups_Index[[2]] = cmap[,c(2,3)][cmap[,c(2,3)]>0]
	parGroups_Index[[3]] = cmap[,c(5,6)][cmap[,c(5,6)]>0]
	parGroups_Index[[4]] = cmap[,8][cmap[,8]>0][1:7]
	parGroups_Index[[5]] = cmap[,12:13][cmap[,12:13]>0]
	parGroups_Index[[6]] = unique(cmap[,14:20][cmap[,14:20]>0])
	parGroups_Index[[7]] = cmap[,21:31][cmap[,21:31]>0]
	#----------------------------------------------------------------------------------------------------------- 


	parChains = list()  # This list stores the MCMC parameter vector chain for each parameter group.  
	for(l in parGroups){ parChains[[l]] = matrix(0,steps,length(parGroups_Index[[l]])) }
	
	Sigma = list()  # This list stores the proposal covariance matrices for each parameter group.
	for(l in parGroups){ Sigma[[l]] = diag( length(parGroups_Index[[l]]) ) }

	tau = rep(10^(-5),7)  # These are used to tune the proposals during the burnin period.
	prev_logPosteriorDens = -Inf  # Initialize.
	C = -par[splineParam[1]]
	D = par[splineParam[2]] + C
	

	### The algorithm... ----------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------
	accept = rep(0,length(parGroups))
	for(ttt in 1:steps){

		### Update each parameter group, in turn.
		for(l in parGroups){ 
			
			# The indices of 'par' corresponding to transition 'l'.
			tempInd = parGroups_Index[[l]]
			
			
			# Proposal scheme -----------------------------------------------------------------------------------
			if(ttt < burnin){
				
				# During the burnin period, only propose parameters which satisfy the constraints.
				constraints_satisfied = FALSE
				count = 1
				while(constraints_satisfied == FALSE & count < 10000){ 
					hold = ProposalFun(par, tempInd, C, D, tau, Sigma, l, designMat, cubicSplineMat, splineMat)
					proposal = hold[[1]]
					constraints_satisfied = hold[[2]]
					
					# Expedite the process of proposing a valid parameter vector.
					if(count %% 100 == 0){  tau[l] = (.5^2)*tau[l]  }
					count = count + 1
				}
			} else{
				
				hold = ProposalFun(par, tempInd, C, D, tau, Sigma, l, designMat, cubicSplineMat, splineMat)
				proposal = hold[[1]]
				constraints_satisfied = hold[[2]]
			}
			#----------------------------------------------------------------------------------------------------
			

			# If the parameter constraits are satisfied, then compute the log Metropolis-Hastings ratio ---------
			if(constraints_satisfied == TRUE){
				
				# Updated spline parameters are only a1 and a2. NO PRIORS ON THE REMAINING CONTROL POINTS.
				logPosteriorDens = fn(proposal) + dmvnorm(proposal[-splineParam[3:11]], 
														 prior.means[-splineParam[3:11]], 
														 diag(prior.std[-splineParam[3:11]]^2), log=TRUE)
												 
				logMH_Ratio = logPosteriorDens - prev_logPosteriorDens
				
			} else{  logMH_Ratio = -Inf  } 
			#----------------------------------------------------------------------------------------------------
			
			
			# Accept or reject the proposed parameters ----------------------------------------------------------
			if( log(runif(1,0,1)) < logMH_Ratio ){ 
				
				# The updates if the proposal is accepted.
				prev_logPosteriorDens = logPosteriorDens
				par = proposal
				accept[l] = accept[l] + 1 
		
			} else{  
				C = -par[splineParam[1]]     # a1.
				D = par[splineParam[2]] + C  # a2. 
			}
			
			# Record the MCMC step.
			parChains[[l]][ttt,] = par[tempInd]
			#----------------------------------------------------------------------------------------------------
			
			
			# Proposal tuning scheme ----------------------------------------------------------------------------
			if(ttt < burnin){
				# During the burnin period, update the proposal covariance in each step to capture the 
				# relationships within the parameters vectors for each transition.  This helps with mixing.
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
		
		# If 10,000 sequential proposals failed to satisfy the constraints, then it's likely the initial values. 
		if(count == 10000){  stop('The initial parameter values do not satisfy the constraints.')  }
	}
	#------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------------
	
	
	accept = accept / (steps-burnin)
	
	list(loglik=prev_logPosteriorDens, tau=tau, accept=accept, par=par, parGroups_Index=parGroups_Index, 
			parChains=parChains, Sigma=Sigma)
}











ProposalFun <- function(par, tempInd, C, D, tau, Sigma, l, designMat, cubicSplineMat, splineMat){
	
	proposal = par  # Initialize.
	
	# Propose parameter(s) as in random walk MH 
	# The parameter group with the spline coefficients requires a special update.
	if(l==4){
		rawProposal = par[tempInd] + mvrnorm(1,rep(0,length(tempInd)),tau[l]*Sigma[[l]])
		C = -rawProposal[6] # a1.
		D = rawProposal[7] + C # a2.

		aVec = seq(-C,by=D,length=11) # Currently there are 11 control points for the b-spline.
		aVec[aVec>0] = 0
		
		proposal[c(tempInd[1:5],splineParam)] = c(rawProposal[1:5], aVec)

	} else{ proposal[tempInd] = par[tempInd] + mvrnorm(1,rep(0,length(tempInd)),tau[l]*Sigma[[l]]) }
	
	
	# Constrain the age coefficients on the transition rates to be nonnegative.
	positiveAgeCoefs = ( prod( proposal[c(13,18,23,28,33,38,43,59,64)] >= 0 ) == 1 )

	# Constrain the control points for the cubic spline on age, for state 1 to state 2, to be non-decreasing.
	cubicSplineConstr = (  proposal[4]<=proposal[5] & proposal[5]<=proposal[6] & proposal[6]<=proposal[7]
					 	 & proposal[7]<=proposal[8] & proposal[8]<=proposal[9] & proposal[9]<=proposal[10]
						 & proposal[10]<=proposal[11] )
	
	# Constrain the MMSE intercepts for states 1-6.
	MMSEconstr = (  proposal[68]>=proposal[74] & proposal[68]>=proposal[75] & proposal[74]>=proposal[76]
			 	  & proposal[75]>=proposal[76] & proposal[76]>=proposal[77] & proposal[77]>=proposal[78] )
	
	# Constrain the amyloid and neurodegeneration response means to be properly separated.
	pibMeansConstr = ( proposal[80] < proposal[81] ) 
	thicknessMeansConstr = ( proposal[83] > proposal[84] )
	
	# Transition rate from A-N+ to A+N+ should be at least as large as the rate from A-N- to A+N-. 
	S1_to_S2 = c( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11)
	S3_to_S4 = c(22,23,24,25,26)
	constr1 = ( prod(designMat %*% proposal[S3_to_S4] >= cubicSplineMat %*% proposal[S1_to_S2]) == 1 )

	# Transition rate from DemA- to DemA+ should be at least as large as the rate from A-N+ to A+N+.
	S5_to_S6 = c(37,38,39,40,41)
	constr2 = ( prod(designMat %*% proposal[S5_to_S6] >= designMat %*% proposal[S3_to_S4]) == 1 )
	
	# Transition rate from A+N- to A+N+ should be at least as large as the rate from A-N- to A-N+.  
	S2_to_S4 = c(17,18,19,20,21)
	S1_to_S3 = c(12,13,14,15,16)
	constr3 = ( prod(designMat %*% proposal[S2_to_S4] >= designMat %*% proposal[S1_to_S3]) == 1 )
	
	# Transition rate from A+N+ to DemA+ should be at least as large as the rate from A-N+ to DemA-.
	S4_to_S6 = c(32,33,34,35,36)
	S3_to_S5 = c(27,28,29,30,31)
	constr4 = ( prod(designMat %*% proposal[S4_to_S6] >= designMat %*% proposal[S3_to_S5]) == 1 )
	
	# Transition rate from DemA- to Dead should be at least as large as the rate from NonDem to Dead.
	S5_to_S7 = c(58,59,60,61,62)
	S1234_to_S7 = c(42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57)
	constr5 = ( prod(designMat %*% proposal[S5_to_S7] >= splineMat %*% proposal[S1234_to_S7]) == 1 )
	
	# Transition rate from DemA+ to Dead should be at least as large as the rate from DemA- to Dead.
	S6_to_S7 = c(63,64,65,66,67)
	constr6 = ( prod(designMat %*% proposal[S6_to_S7] >= designMat %*% proposal[S5_to_S7]) == 1 )
	
	constraints_satisfied <- ( constr1 & constr2 & constr3 & constr4 & constr5 & constr6 & C>0 & D>0 & C>D
							   & MMSEconstr & positiveAgeCoefs & pibMeansConstr & thicknessMeansConstr 
							   & cubicSplineConstr )
	   
	return( list( proposal, constraints_satisfied) )
							  
}




