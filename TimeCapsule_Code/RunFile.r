library(hmm)
source('SetupFunc.r')

args <- commandArgs(TRUE)
set.seed(as.integer(args[1]))
sensitivity_test <- args[4]
NCORES <- as.integer(args[5])
Dir <- args[6]

if(args[2]=='Simulation'){
	load(paste0(Dir,'demData',args[1],'.rda'))
	DATA <- demData
	STEPS  <- 15000
	BURNIN <- 10000
	
} else if(args[2]=='RealData'){	
	load('synthetic_MCSA_data.rda')
	DATA <- useData
	STEPS  <- 30000
	BURNIN <- 10000
	
}



states <- c("A-N-", "A+N-", "A-N+", "A+N+", "DemA-", "DemA+", "Dead")

qmat <- matrix(c(0, .1, .1,  0,  0,  0, .09,
                 0,  0,  0, .1,  0,  0, .09,
                 0,  0,  0, .1, .1,  0, .09,
                 0,  0,  0,  0,  0, .1, .09,
                 0,  0,  0,  0,  0, .1,  .1,
				 0,  0,  0,  0,  0,  0,  .1,
                 0,  0,  0,  0,  0,  0,   0), ncol=7, byrow=TRUE)
dimnames(qmat) <- list(states, states)


qcoef <- data.frame(state1 = c(  rep(1,12),  rep(1,5),  rep(2,5),  rep(3,5), rep(3,5), rep(4,5), rep(5,5),
								 rep(1,16), rep(2,16), rep(3,16), rep(4,16), rep(5,5), rep(6,5)  ),
								 
				    state2 = c(  rep(2,12),  rep(3,5),  rep(4,5),  rep(4,5), rep(5,5), rep(6,5), rep(6,5), 
								 rep(7,16), rep(7,16), rep(7,16), rep(7,16), rep(7,5), rep(7,5)  ),
								 
                    term   = c('(Intercept)','male','educ','apoe4','c1_amyl','c2_amyl','c3_amyl','c4_amyl',
													 'c5_amyl','c6_amyl','c7_amyl','c8_amyl',
							   '(Intercept)','iage','male','educ','apoe4',
							   '(Intercept)','iage','male','educ','apoe4',
							   '(Intercept)','iage','male','educ','apoe4', 
							   '(Intercept)','iage','male','educ','apoe4',
							   '(Intercept)','iage','male','educ','apoe4',
							   '(Intercept)','iage','male','educ','apoe4', 
							   '(Intercept)','iage','male','educ','apoe4','c1','c2','c3','c4','c5','c6','c7',
							   											  'c8','c9','c10','c11',
							   '(Intercept)','iage','male','educ','apoe4','c1','c2','c3','c4','c5','c6','c7',
							   											  'c8','c9','c10','c11', 
							   '(Intercept)','iage','male','educ','apoe4','c1','c2','c3','c4','c5','c6','c7',
							   											  'c8','c9','c10','c11',
							   '(Intercept)','iage','male','educ','apoe4','c1','c2','c3','c4','c5','c6','c7',
							   											  'c8','c9','c10','c11', 
							   '(Intercept)','iage','male','educ','apoe4',
							   '(Intercept)','iage','male','educ','apoe4'),
                    coef   = c( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,
							   12,13,14,15,16,
							   17,18,19,20,21,
							   22,23,24,25,26,
							   27,28,29,30,31,
							   32,33,34,35,36,
							   37,38,39,40,41,
							   42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,
							   42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,
							   42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,
							   42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,
							   58,59,60,61,62,
							   63,64,65,66,67),

					init   = c( 0.000000000,-0.331664059, 0.014273922, 0.945136344,-6.580919931,-6.501250885,
															     				   -6.427228628,-4.673031945,
																			       -3.965086053,-3.111756265,
																			       -1.950950409, 3.588549810,
						       -5.000000000, 0.094236594, 0.141479641,-0.095807126, 0.486230191,
						       -4.200000000, 0.093838248,-0.116771101,-0.098580487, 0.399355964,
						       -3.881605509, 0.065857892,-0.330382263, 0.019152541, 0.444423684,
						       -4.148148694, 0.075589995, 0.248030366, 0.033380640, 0.765503215,
						       -2.473865801, 0.045208611, 0.082975736, 0.054272656, 0.256422760,
						       -2.739931535, 0.047242201, 0.050003054, 0.033386975, 0.239889316,
						       -4.329458534, 0.110780340, 0.387452143,-0.042883267,-0.008324821,
							   			    -1.242275027,-1.112106935,-0.981938842,-0.851770750,-0.721602657,
							 			    -0.591434565,-0.461266472,-0.331098380,-0.200930287,-0.070762195,
											 0.000000000,
  						       -4.329458534, 0.110780340, 0.387452143,-0.042883267,-0.008324821,
  							   			    -1.242275027,-1.112106935,-0.981938842,-0.851770750,-0.721602657,
  							 			    -0.591434565,-0.461266472,-0.331098380,-0.200930287,-0.070762195,
  											 0.000000000,
  						       -4.329458534, 0.110780340, 0.387452143,-0.042883267,-0.008324821,
  							   			    -1.242275027,-1.112106935,-0.981938842,-0.851770750,-0.721602657,
  							 			    -0.591434565,-0.461266472,-0.331098380,-0.200930287,-0.070762195,
  											 0.000000000,
  						       -4.329458534, 0.110780340, 0.387452143,-0.042883267,-0.008324821,
  							   			    -1.242275027,-1.112106935,-0.981938842,-0.851770750,-0.721602657,
  							 			    -0.591434565,-0.461266472,-0.331098380,-0.200930287,-0.070762195,
  											 0.000000000,
						       -3.500000000, 0.120000000, 0.360000000,-0.007003006, 0.104718992,
						       -1.100000000, 0.060000000, 0.224925174, 0.003157103, 0.236539477), 
							   stringsAsFactors=FALSE )



rcoef <- data.frame(response = c(rep(1,37),rep(2,3),rep(3,3),4,4),
					lp   = c(1,1,1,1,1,1, 
							 2,2,2,2,2,2,
							 3,3,3,3,3,3,
							 4,4,4,4,4,4,
							 5,5,5,5,5,5,
							 6,6,6,6,6,6,
							 7,
							 
							 8,9,10,
							 
							 11,12,13,
							 
							 14,15),
					
					term = c('(Intercept)','age','male','educ','apoe4','ntests', # State 1 (Intercept).
							 '(Intercept)','age','male','educ','apoe4','ntests', # State 2 (Intercept).
							 '(Intercept)','age','male','educ','apoe4','ntests', # State 3 (Intercept).
							 '(Intercept)','age','male','educ','apoe4','ntests', # State 4 (Intercept).
							 '(Intercept)','age','male','educ','apoe4','ntests', # State 5 (Intercept).
							 '(Intercept)','age','male','educ','apoe4','ntests', # State 6 (Intercept).
							 '(Intercept)',
							 
							 rep('(Intercept)',8)),
							 
					coef = c( 1,2,3,4,5,6, 
							  7,2,3,4,5,6,
							  8,2,3,4,5,6,
							  9,2,3,4,5,6,
							 10,2,3,4,5,6,
							 11,2,3,4,5,6,
							 12,
							 
							 13,14,15,
							 
							 16,17,18,
							 
							 19,20),
							
					init = c( 0.538414931,-0.014554096,-0.089207557, 0.058688151,-0.004627584, 0.057826362,
  						      0.178435764,-0.014554096,-0.089207557, 0.058688151,-0.004627584, 0.057826362,
  							 -0.672421219,-0.014554096,-0.089207557, 0.058688151,-0.004627584, 0.057826362,
  							 -1.827065032,-0.014554096,-0.089207557, 0.058688151,-0.004627584, 0.057826362,
  							 -2.341332296,-0.014554096,-0.089207557, 0.058688151,-0.004627584, 0.057826362,
							 -5.208797471,-0.014554096,-0.089207557, 0.058688151,-0.004627584, 0.057826362,
  							 -1.261838927,
							 
							 -1.084010116, 0.110684536,-2.240498159,
							 
							  2.655601472, 2.498948044,-3.869343130,
							  
							 -4.712799832,-2.163822374), stringsAsFactors=FALSE )
							
							

pcoef <- data.frame(lp   = c(1,2,3),		 
					term = c('(Intercept)','(Intercept)','(Intercept)'),
					coef = c(1,2,3),
					init = c(-3.9835,-4.9499,-6.5185), stringsAsFactors=FALSE )
							 
							 

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



lpibM <-  c( -1.3, -.5, log((.4/3)^2))
lpibSd <- c(   .2,  .2,             2)

thickM <-  c( 3.14, 2.34, log((.4/3)^2))
thickSd <- c(  .2,    .2,             2)

# Death rate adjustment for women.
# 507.1/590.4 = 0.859
# log(exp(-4.26)* 0.859) = -4.411986

# Death rate adjustment for men.
# 724.8/817.3 = 0.887
# log(exp(-4.26+.44)* 0.887) = -3.93991
# -3.93991--4.411986 = 0.472076

means <- matrix(c( NA,-3,-3,-3,-3,-3,-3,-4.41,-4,-4,   rep(-.28,4),-7.3,-7.3,-.7,  lpibM,thickM,-3,-3,-3.5,-6,-6,
		 		    0,.1,.1,.1,.1,.1,.1, .094,.1,.1,    0,NA,NA,NA,  NA,  NA, NA,  rep(NA,11),
				    0, 0, 0, 0, 0, 0, 0,  .47, 0, 0,    0,NA,NA,NA,  NA,  NA, NA,  rep(NA,11),
				    0, 0, 0, 0, 0, 0, 0,    0, 0, 0,    0,NA,NA,NA,  NA,  NA, NA,  rep(NA,11),
				   -5, 0, 0, 0, 0, 0, 0,    0, 0, 0,    0,NA,NA,NA,  NA,  NA, NA,  rep(NA,11),
				   -4,NA,NA,NA,NA,NA,NA, -.75,NA,NA,    0,NA,NA,NA,  NA,  NA, NA,  rep(NA,11),
				   -3,NA,NA,NA,NA,NA,NA, -.60,NA,NA,   NA,NA,NA,NA,  NA,  NA, NA,  rep(NA,11),
				   -2,NA,NA,NA,NA,NA,NA,    0,NA,NA,   NA,NA,NA,NA,  NA,  NA, NA,  rep(NA,11),
				   -1,NA,NA,NA,NA,NA,NA,    0,NA,NA,   NA,NA,NA,NA,  NA,  NA, NA,  rep(NA,11),
				    0,NA,NA,NA,NA,NA,NA,    0,NA,NA,   NA,NA,NA,NA,  NA,  NA, NA,  rep(NA,11),
				    1,NA,NA,NA,NA,NA,NA,    0,NA,NA,   NA,NA,NA,NA,  NA,  NA, NA,  rep(NA,11),
				    2,NA,NA,NA,NA,NA,NA,    0,NA,NA,   NA,NA,NA,NA,  NA,  NA, NA,  rep(NA,11),
				   NA,NA,NA,NA,NA,NA,NA,    0,NA,NA,   NA,NA,NA,NA,  NA,  NA, NA,  rep(NA,11),
				   NA,NA,NA,NA,NA,NA,NA,    0,NA,NA,   NA,NA,NA,NA,  NA,  NA, NA,  rep(NA,11),
				   NA,NA,NA,NA,NA,NA,NA,    0,NA,NA,   NA,NA,NA,NA,  NA,  NA, NA,  rep(NA,11),
				   NA,NA,NA,NA,NA,NA,NA,    0,NA,NA,   NA,NA,NA,NA,  NA,  NA, NA,  rep(NA,11)
				   ), 16, 28, byrow=TRUE)
prior.means <- c(means[!is.na(means)])		
 
if(sensitivity_test=='no_scale'){								 
stdev <- matrix(c( NA,  1,  1,  1,  1,  1,  1,  .1,  1,  1,    rep(.75,4), 3, 3, 2,  lpibSd,thickSd,1,1,.25,1,1,
				    1,.05,.05,.05,.05,.05,.05, .01,.05,.05,    1,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				   .1,  1,  1,  1,  1,  1,  1, .05,  1,  1,    1,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				    1, .1, .1, .1, .1, .1, .1,  .1, .1, .1,    1,NA,NA,NA,NA,NA,NA,  rep(NA,11),
			 	    1,  1,  1,  1,  1,  1,  1,   1,  1,  1,    1,NA,NA,NA,NA,NA,NA,  rep(NA,11),
			 	    2, NA, NA, NA, NA, NA, NA,.375, NA, NA,    1,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				    2, NA, NA, NA, NA, NA, NA,  .3, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				    2, NA, NA, NA, NA, NA, NA,   0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				    2, NA, NA, NA, NA, NA, NA,   0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				    3, NA, NA, NA, NA, NA, NA,   0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				    3, NA, NA, NA, NA, NA, NA,   0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				    3, NA, NA, NA, NA, NA, NA,   0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				   NA, NA, NA, NA, NA, NA, NA,   0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				   NA, NA, NA, NA, NA, NA, NA,   0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				   NA, NA, NA, NA, NA, NA, NA,   0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				   NA, NA, NA, NA, NA, NA, NA,   0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11)
				   ), 16, 28, byrow=TRUE)
} else if(sensitivity_test=='10x_scale'){			   
stdev <- matrix(c( NA,  1,  1,  1,  1,  1,  1,  .01,  1,  1,    rep(.75,4), 3, 3, 2,  lpibSd,thickSd,1,1,.25,1,1,
				    1,.05,.05,.05,.05,.05,.05, .001,.05,.05,    1,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				   .1,  1,  1,  1,  1,  1,  1, .005,  1,  1,    1,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				    1, .1, .1, .1, .1, .1, .1,  .01, .1, .1,    1,NA,NA,NA,NA,NA,NA,  rep(NA,11),
			 	    1,  1,  1,  1,  1,  1,  1,   .1,  1,  1,    1,NA,NA,NA,NA,NA,NA,  rep(NA,11),
			 	    2, NA, NA, NA, NA, NA, NA,.0375, NA, NA,    1,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				    2, NA, NA, NA, NA, NA, NA,  .03, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				    2, NA, NA, NA, NA, NA, NA,    0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				    2, NA, NA, NA, NA, NA, NA,    0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				    3, NA, NA, NA, NA, NA, NA,    0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				    3, NA, NA, NA, NA, NA, NA,    0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				    3, NA, NA, NA, NA, NA, NA,    0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				   NA, NA, NA, NA, NA, NA, NA,    0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				   NA, NA, NA, NA, NA, NA, NA,    0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				   NA, NA, NA, NA, NA, NA, NA,    0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11),
				   NA, NA, NA, NA, NA, NA, NA,    0, NA, NA,   NA,NA,NA,NA,NA,NA,NA,  rep(NA,11)
			   ), 16, 28, byrow=TRUE) * 10
}
prior.std <- c(stdev[!is.na(stdev)])




if(args[3] == 'Bayesian'){
	
	source('mcmcRoutine.r')
	
	splineParam = 47:57
	designMat <- as.matrix( cbind( '(Intercept)'=rep(1,nrow(DATA)), DATA[,c('iage','male','educ','apoe4')]) )
	cubicSplineMat <- as.matrix( DATA[,c('male','educ','apoe4','c1_amyl','c2_amyl','c3_amyl','c4_amyl','c5_amyl',
	                                     'c6_amyl','c7_amyl','c8_amyl')])
	splineMat <- as.matrix( cbind( '(Intercept)'=rep(1,nrow(DATA)), DATA[,c('iage','male','educ','apoe4','c1',
																			'c2','c3','c4','c5','c6','c7','c8',
																			'c9','c10','c11')]) )
	
	Output <- hmm(hbind(age, mmse, lpib, thickness, DemStatus) ~ 1 + iage + age + male + educ + apoe4 + ntests + 
					c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10 + c11 + c1_amyl + c2_amyl + c3_amyl + 
					c4_amyl + c5_amyl + c6_amyl + c7_amyl + c8_amyl, data=DATA, scale=FALSE,
					mc.cores=NCORES, entry=c(1,1,1,1,0,0,0), id=ptnum, qmatrix=qmat, qcoef=qcoef, death=7, 
					rfun=list(MMSEResp, lpibResp, thickResp, demResp), otype=obstype, rcoef=rcoef, pfun=pfun, 
					pcoef=pcoef, mfun=mcmcRoutine, mpar=list(cmap=cmap, splineParam=splineParam, steps=STEPS, 
					burnin=BURNIN, prior.means=prior.means, prior.std=prior.std, designMat=designMat,
					splineMat=splineMat, cubicSplineMat=cubicSplineMat))

} else if(args[3] == 'ML'){
	
	Output <- hmm(hbind(age, mmse, lpib, thickness, DemStatus) ~ 1 + iage + age + male + educ + apoe4 + ntests + 
					c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10 + c11 + c1_amyl + c2_amyl + c3_amyl + 
					c4_amyl + c5_amyl + c6_amyl + c7_amyl + c8_amyl, data=DATA, scale=FALSE,
					mc.cores=NCORES, entry=c(1,1,1,1,10^-8,10^-8,10^-8), id=ptnum, qmatrix=qmat, qcoef=qcoef, 
					death=7, rfun=list(MMSEResp, lpibResp, thickResp, demResp), otype=obstype, rcoef=rcoef, 
					pfun=pfun, pcoef=pcoef, mfun=optim, mpar=list(control=list(fnscale= -1, maxit=100), 
					method='BFGS', gr= "hmmgrad", hessian=TRUE))
}

save(Output, file=paste0(Dir,'Output',args[1],'.rda'))

