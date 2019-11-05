# This file only containts stand-alone R functions for implementing the HMM on the CAV simuluated data.
# These functions are required for 'RunFile_cav.r'.






errfun2 <- function(y, nstate, eta, gradient) {
    # true states 1, 2, and 3 have separate linear predictors
    temp1 <- hmulti(y, nstate, eta[,1], gradient,
                    statemap= rbind(1:2, 0,0,0))
    temp2 <- hmulti(y, nstate, eta[,2:3], gradient,
                    statemap= rbind(0, c(2,1,3), 0, 0))
    temp3 <- hmulti(y, nstate, eta[,4], gradient,
                    statemap= rbind(0, 0, c(3,2), 0))
    pmat <- rbind(temp1[1,], temp2[2,], temp3[3,], ifelse(y==4,1,0))
    if (gradient) {
        gmat <- array(c(attr(temp1, 'gradient'), attr(temp2, 'gradient'),
                        attr(temp3, 'gradient')), dim=c(nstate, length(y), 4))
        attr(pmat, "gradient") <- gmat
    }
    pmat
}





misclassResp <- function( y, nstate, eta, gradient){  #, map=c(1,2,3,4), setup=misclass_setup) {
	
	if(gradient){
		errfun2(y, nstate, eta, gradient)
	
	} else{
		
		# The misclassification matrix is 'true state' by 'observed state'
		temp = matrix(c(            1, exp(eta[1,1]),             0, 0,
		       		    exp(eta[1,2]),             1, exp(eta[1,3]), 0,
				 	                0, exp(eta[1,4]),             1, 0,
					                0,             0,             0, 1), ncol=4, byrow=TRUE)
		temp = temp / rowSums(temp)

		# nstate by n matrix yielding P(observed state | true state = state j), where j in {1,...,nstate}.
		emat <- temp[,y, drop=FALSE]
		rownames(emat) = c('Well','Mild','Severe','Death')
	
		return( emat ) 
	}
    
}





pfun <- function( nstate, eta, gradient) {
	
	if(gradient){
		hmminit(nstate, eta, gradient)
	
	} else{
		initProbs = c(1, exp(eta[1]), exp(eta[2]), 0)
		return( initProbs / sum(initProbs) )
	}
}