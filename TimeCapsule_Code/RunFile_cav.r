library(hmm)
library(msm)
source('mcmcRoutine_cav.r')
source('SetupFunc_cav.r')

args <- commandArgs(TRUE)
set.seed(as.integer(args[1]))
NCORES <- as.integer(args[2])
Dir <- args[3]
load(paste0(Dir,'cavData',args[1],'.rda'))




qmat <- matrix(c( 0,exp(-2),      0,exp(-2),
                  0,      0,exp(-2),exp(-2),
                  0,      0,      0,exp(-2),
                  0,      0,      0,      0), ncol=4, byrow=TRUE)
dimnames(qmat) <- list( c('Well', 'Mild','Severe','Death'), c('Well', 'Mild','Severe','Death'))


qcoef <- data.frame(state1 = c(  rep(1,3),  rep(1,3),  rep(2,3), rep(2,3), rep(3,3)),
				    state2 = c(  rep(2,3),  rep(4,3),  rep(3,3), rep(4,3), rep(4,3)),
                    term   = c('(Intercept)','iyears','sex', 
							   '(Intercept)','iyears','sex',
							   '(Intercept)','iyears','sex',
							   '(Intercept)','iyears','sex',
							   '(Intercept)','iyears','sex'),
                    coef   = c( 1, 2, 3,
							    4, 5, 6,
							    7, 8, 9,
							   10,11,12,
							   13,14,15),
					init   = c(-2, 0, 0,
							   -2, 0, 0,
							   -2, 0, 0,
							   -2, 0, 0,
							   -2, 0, 0), stringsAsFactors=FALSE) 


rcoef <- data.frame(response = c(1,1,1,1),
					lp	   	 = c(1,2,3,4),
					term     = c('(Intercept)','(Intercept)','(Intercept)','(Intercept)'),
					coef 	 = c(1,2,3,4),
					init 	 = rep(-3,4), stringsAsFactors=FALSE )  


pcoef <- data.frame(lp   = c(1,2),		 
					term = c('(Intercept)','(Intercept)'),
					coef = c(1,2),
					init = c(-4.5,-5), stringsAsFactors=FALSE )


cmap <- matrix(c( 1, 4, 7,10,13,  16,17,18,19,20,21,
 				  2, 5, 8,11,14,   0, 0, 0, 0, 0, 0,
 				  3, 6, 9,12,15,   0, 0, 0, 0, 0, 0), 3, 11, byrow=TRUE)






#----------------------------------------------------------------------------------------------------------------
# Run the Bayesian implementation -------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

means <- matrix(c(-2,-2,-2,-2,-2,  -3,-3,-3,-3,-4.5,-5,
				   0, 0, 0, 0, 0,  NA,NA,NA,NA,  NA,NA,
				   0, 0, 0, 0, 0,  NA,NA,NA,NA,  NA,NA), 3, 11, byrow=TRUE)
prior.means <- c(means[!is.na(means)])		

			 
#stdev <- matrix(c( 1, 1, 1, 1, 1,  1.5,1.5,1.5,1.5,1.5,1.5,
#				   .5,.5,.5,.5,.5,   NA, NA, NA, NA, NA, NA,
#				    1, 1, 1, 1, 1,   NA, NA, NA, NA, NA, NA), 3, 11, byrow=TRUE)
stdev <- matrix(c( 2, 2, 2, 2, 2,    3,  3,  3,  3,  3,  3,
				   1, 1, 1, 1, 1,   NA, NA, NA, NA, NA, NA,
				   2, 2, 2, 2, 2,   NA, NA, NA, NA, NA, NA), 3, 11, byrow=TRUE)
prior.std <- c(stdev[!is.na(stdev)])


Output_hmm <- hmm(hbind(years, state) ~ 1 + iyears + sex, data=cavData, mc.cores=NCORES, entry=c(1,1,1,0), 
				  id=ptnum, rfun=list(misclassResp), rcoef=rcoef, pfun=pfun, pcoef=pcoef, qmatrix=qmat, 
				  qcoef=qcoef, death=4, otype=obstype_hmm, mfun=mcmcRoutine, scale=FALSE,
				  mpar=list(cmap=cmap, prior.means=prior.means, prior.std=prior.std, steps=15000, burnin=5000))

save(Output_hmm, file=paste0(Dir,'Output_hmm',args[1],'.rda'))

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
# Run the maximum likelihood implementation ---------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

Output_mle <- hmm(hbind(years, state) ~ 1 + iyears + sex, data=cavData, mc.cores=NCORES, entry=c(1,1,1,0), 
                  death=4, id=ptnum, rfun=list(misclassResp), rcoef=rcoef, pfun=pfun, pcoef=pcoef, qmatrix=qmat, 
				  qcoef=qcoef, otype=obstype_hmm, scale=FALSE, mfun=hmmscore, 
				  mpar=list(gr="hmmboth", iter=10000))
				  
save(Output_mle, file=paste0(Dir,'Output_mle',args[1],'.rda'))

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
# Run the msm implementation ------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

emat = matrix(c(      1, exp(-3),       0, 0,
       		    exp(-3),       1, exp(-3), 0,
		 	          0, exp(-3),       1, 0,
			          0,       0,       0, 1), ncol=4, byrow=TRUE)
emat = emat / rowSums(emat)
dimnames(emat) <- list( c('Well','Mild','Severe','Death'), c('Well','Mild','Severe','Death'))


Output_msm <- msm(state ~ years, subject=ptnum, data=cavData, qmatrix=qmat, covariates= ~ 1 + iyears + sex, 
				  center=FALSE, covinits=list(iyears=c(0,0,0,0,0),sex=c(0,0,0,0,0)), obstrue=obstrue, 
				  ematrix=emat, initprobs=c(1, exp(-4.5), exp(-5), 0), est.initprobs=TRUE, deathexact=4, 
				  censor=99, censor.states=1:3, method='BFGS', control=list(fnscale=4000, maxit=10000))   

save(Output_msm,file=paste0(Dir,'Output_msm',args[1],'.rda'))

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------



# Generate the true parameter values from the real data set.
#load('useData_SimSet1.rda')
#Output_msm <- msm(state ~ years, subject=PTNUM, data=useData_cav, qmatrix=qmat, obstrue=obstrue, ematrix=emat, 
#			  covariates= ~ 1 + iyears + sex, center=FALSE, deathexact=4, censor=99, censor.states=1:3, 
#			  method='BFGS', control=list(fnscale=4000, maxit=10000))


