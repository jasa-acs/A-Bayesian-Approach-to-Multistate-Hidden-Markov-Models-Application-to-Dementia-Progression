# This file only containts stand-alone R functions for implementing the HMM on the MCSA and synthetic MCSA data.
# These functions are required for 'RunFile.r'.



MMSEResp <- function(y,nstate,eta,gradient){
	
	# log(sigma^2) is the parameter for inference because a Gaussian prior is meaningful, and a Gaussian proposal is appropriate.  Recall that all parameters have Gaussain priors and updates.
	sigma_sq = exp(eta[1,7])
	
	row1 = exp( -.5* ( y - eta[,1] )^2 / sigma_sq ) / sqrt( 2*pi*sigma_sq )
	row1[is.na(row1)] = 1
	
	row2 = exp( -.5* ( y - eta[,2] )^2 / sigma_sq ) / sqrt( 2*pi*sigma_sq )
	row2[is.na(row2)] = 1
	
	row3 = exp( -.5* ( y - eta[,3] )^2 / sigma_sq ) / sqrt( 2*pi*sigma_sq )
	row3[is.na(row3)] = 1
	
	row4 = exp( -.5* ( y - eta[,4] )^2 / sigma_sq ) / sqrt( 2*pi*sigma_sq )
	row4[is.na(row4)] = 1
	
	row5 = exp( -.5* ( y - eta[,5] )^2 / sigma_sq ) / sqrt( 2*pi*sigma_sq )
	row5[is.na(row5)] = 1
	
	row6 = exp( -.5* ( y - eta[,6] )^2 / sigma_sq ) / sqrt( 2*pi*sigma_sq )
	row6[is.na(row6)] = 1
	
	# nstate by n matrix yielding P(y_i | s_i = state j), where j in {1,...,nstate}.
	rmat = as.matrix( rbind(row1, 
							row2, 
							row3, 
							row4, 
							row5, 
							row6, 
							rep(0,length(y))) )
	rownames(rmat) = c("A-N-", "A+N-", "A-N+", "A+N+", "DemA-", "DemA+", "Dead")
	
	
	# If using a gradient based procedure, then assign gradient attributes --------------------------------------
	if(gradient){
		gmat = array(0., dim=c(nstate, length(y), ncol(eta)))
		
		# Update each of the 7 rows of the 7 gradient arrays.  However, the 7th rows are all zeros so k = 1:6.  Note that the only nonzero row for the jth partial derivative is the jth row, for j \in {1,...,6}, and each of the rows for the 7th partial derivative have the same form.  That is why there is only one for-loop.
		for(k in 1:6){
			
			x = (y - eta[,k])
			gmat[k,,k] = x * exp( - .5* log(sigma_sq) - .5* x^2 / sigma_sq ) / (sigma_sq*sqrt(2*pi))
			gmat[k,,7] = exp( - .5* log(sigma_sq) - .5* x^2 / sigma_sq ) *( .5 * x^2 / sigma_sq - .5 )/sqrt(2*pi)
		}
		
		attr( rmat, 'gradient') = gmat
	}
	#------------------------------------------------------------------------------------------------------------
	
	
	return( rmat ) 
}





lpibResp <- function(y,nstate,eta,gradient) {
	
	# log(sigma^2) is the parameter for inference because a Gaussian prior is meaningful, and a Gaussian proposal is appropriate.  Recall that all parameters have Gaussain priors and updates.
	sigma_sq = exp(eta[1,3])
	
	A_minus = exp( -( y - eta[,1] )^2 / (2*sigma_sq) ) / sqrt( 2*pi*sigma_sq )
	
	A_plus = exp( -( y - eta[,2] )^2 / (2*sigma_sq) ) / sqrt( 2*pi*sigma_sq )
	
	# nstate by n matrix yielding P(y_i | s_i = state j), where j in {1,...,nstate}.
 	rmat = as.matrix( rbind(A_minus, 
 							A_plus, 
 							A_minus, 
 							A_plus, 
 							A_minus, 
 							A_plus, 
 							rep(0,length(y))) )
 	rownames(rmat) = c("A-N-", "A+N-", "A-N+", "A+N+", "DemA-", "DemA+", "Dead")
	
	
	# If using a gradient based procedure, then assign gradient attributes --------------------------------------
	if(gradient){
		gmat = array(0., dim=c(nstate, length(y), ncol(eta)))
		
		# Update each of the 7 rows of the 3 gradient arrays.  However, the 7th rows are all zeros so k = 1:6.
		for(k in 1:6){
			
			if( k %in% c(1,3,5) ){  
				x = (y - eta[,1]) 
				gmat[k,,1] = x * exp( - .5* log(sigma_sq) - .5* x^2 / sigma_sq ) / (sigma_sq*sqrt(2*pi))
				gmat[k,,2] = 0 
			
			} else if( k %in% c(2,4,6) ){  
				x = (y - eta[,2]) 
				gmat[k,,1] = 0
				gmat[k,,2] = x * exp( - .5* log(sigma_sq) - .5* x^2 / sigma_sq ) / (sigma_sq*sqrt(2*pi))
			
			}
			gmat[k,,3] = exp( - .5* log(sigma_sq) - .5* x^2 / sigma_sq ) *( .5 * x^2 / sigma_sq - .5 )/sqrt(2*pi)
		}
		
		attr( rmat, 'gradient') = gmat
	}
	#------------------------------------------------------------------------------------------------------------
	
	
	return( rmat )    
}





thickResp <- function(y,nstate,eta,gradient) {
	
	# log(sigma^2) is the parameter for inference because a Gaussian prior is meaningful, and a Gaussian proposal is appropriate.  Recall that all parameters have Gaussain priors and updates.
	sigma_sq = exp(eta[1,3])
	
	N_minus = exp( -( y - eta[,1] )^2 / (2*sigma_sq) ) / sqrt( 2*pi*sigma_sq )
	
	N_plus = exp( -( y - eta[,2] )^2 / (2*sigma_sq) ) / sqrt( 2*pi*sigma_sq )
 	
	# nstate by n matrix yielding P(y_i | s_i = state j), where j in {1,...,nstate}.
 	rmat = as.matrix( rbind(N_minus, 
 							N_minus, 
 							N_plus, 
 							N_plus, 
 							N_plus, 
 							N_plus, 
 							rep(0,length(y))) )
 	rownames(rmat) = c("A-N-", "A+N-", "A-N+", "A+N+", "DemA-", "DemA+", "Dead")
	
	
	# If using a gradient based procedure, then assign gradient attributes --------------------------------------
	if(gradient){
		gmat = array(0., dim=c(nstate, length(y), ncol(eta)))
		
		# Update each of the 7 rows of the 3 gradient arrays.  However, the 7th rows are all zeros so k = 1:6.
		for(k in 1:6){
			
			if( k %in% c(1,2) ){  
				x = (y - eta[,1]) 
				gmat[k,,1] = x * exp( - .5* log(sigma_sq) - .5* x^2 / sigma_sq ) / (sigma_sq*sqrt(2*pi))
				gmat[k,,2] = 0 
			
			} else if( k %in% c(3,4,5,6) ){  
				x = (y - eta[,2]) 
				gmat[k,,1] = 0
				gmat[k,,2] = x * exp( - .5* log(sigma_sq) - .5* x^2 / sigma_sq ) / (sigma_sq*sqrt(2*pi))
			
			}
			gmat[k,,3] = exp( - .5* log(sigma_sq) - .5* x^2 / sigma_sq ) *( .5 * x^2 / sigma_sq - .5 )/sqrt(2*pi)
		}
		
		attr( rmat, 'gradient') = gmat
	}
	#------------------------------------------------------------------------------------------------------------
	
	
	return( rmat ) 
}





demResp <- function(y,nstate,eta,gradient){
	
	# The misclassification matrix is 'true state' by c('NonDem','Dem').
	temp = matrix(c(             1, exp(eta[1,1]),
				                 1, exp(eta[1,1]), 
					             1, exp(eta[1,1]),
					             1, exp(eta[1,1]),
				     exp(eta[1,2]),             1, 
					 exp(eta[1,2]),             1,
				                 0,             0), ncol=2, byrow=TRUE)
	temp = temp / rowSums(temp) # Rows sum to 1 (R matrices are in column-major order).
	temp[is.na(temp)] <- 0
	colnames(temp) = c('NonDem','Dem')
	rownames(temp) = c('A-N-', 'A+N-', 'A-N+', 'A+N+', 'DemA-', 'DemA+', 'Dead')
	
	rmat = temp[,y, drop=FALSE]
	
	
	# If using a gradient based procedure, then assign gradient attributes --------------------------------------
	if(gradient){
        
		gmat = array(0., dim=c(nstate, length(y), ncol(eta)))
		
		# Partial derivatives of rmat with respect to eta1.
		gmat_eta1 <- matrix(0.,7,2)
		comp1 = - exp(eta[1,1]) * ( 1 + exp(eta[1,1]) )^-2
		comp2 =   exp(eta[1,1]) * ( 1 + exp(eta[1,1]) )^-2
		gmat_eta1[1,] = c( comp1,  comp2)
		gmat_eta1[2,] = c( comp1,  comp2)
		gmat_eta1[3,] = c( comp1,  comp2)
		gmat_eta1[4,] = c( comp1,  comp2)
		#-------------------------------------------------
		
		# Partial derivatives of rmat with respect to eta2.
		gmat_eta2 <- matrix(0.,7,2)
		comp1 =   exp(eta[1,2]) * ( 1 + exp(eta[1,2]) )^-2
		comp2 = - exp(eta[1,2]) * ( 1 + exp(eta[1,2]) )^-2
		gmat_eta2[5,] = c( comp1,  comp2)
		gmat_eta2[6,] = c( comp1,  comp2)
		#-------------------------------------------------
		
		gmat[,,1] = gmat_eta1[,y, drop=FALSE]
		gmat[,,2] = gmat_eta2[,y, drop=FALSE]
	
		attr( rmat, 'gradient') = gmat
	} 
	#------------------------------------------------------------------------------------------------------------
	

	return( rmat ) 
}





pfun <- function(nstate,eta,gradient) {
	
	if(gradient){
		hmminit(nstate, eta, gradient)
	
	} else{
		initProbs = c(1, exp(eta[1]), exp(eta[2]), exp(eta[3]), 0, 0, 0)
		return( initProbs / sum(initProbs) )
	}
}