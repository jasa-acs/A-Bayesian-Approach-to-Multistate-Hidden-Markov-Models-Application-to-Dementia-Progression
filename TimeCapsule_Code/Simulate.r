library(splines)
library(expm)
args <- commandArgs(TRUE)
set.seed(as.integer(args[1]))
Dir <- args[2]


# Choose the discretization of time.
dt <- 1/365


# These are the true parameter values for the uncentered data ( intercept - coef*mean ).
betaMat <- matrix(c(          NA,-0.331664059, 0.014273922, 0.945136344,          NA,
				    -5.000000000, 0.094236594, 0.141479641,-0.095807126, 0.486230191,
				    -4.200000000, 0.093838248,-0.116771101,-0.098580487, 0.399355964,
				    -3.881605509, 0.065857892,-0.330382263, 0.019152541, 0.444423684,
				    -4.148148694, 0.075589995, 0.248030366, 0.033380640, 0.765503215,
				    -2.473865801, 0.045208611, 0.082975736, 0.054272656, 0.256422760,
				    -2.739931535, 0.047242201, 0.050003054, 0.033386975, 0.239889316,
				    -4.329458534, 0.110780340, 0.387452143,-0.042883267,-0.008324821,
				    -4.329458534, 0.110780340, 0.387452143,-0.042883267,-0.008324821,
				    -4.329458534, 0.110780340, 0.387452143,-0.042883267,-0.008324821,
				    -4.329458534, 0.110780340, 0.387452143,-0.042883267,-0.008324821,
				    -3.500000000, 0.120000000, 0.360000000,-0.007003006, 0.104718992,
				    -1.100000000, 0.060000000, 0.224925174, 0.003157103, 0.236539477), ncol=5, byrow=TRUE)
									
betaMat[,1] <- betaMat[,1] - betaMat[,2]*68 - betaMat[,4]*14.03802

initProbs <- c( 1, exp(-3.9835), exp(-4.9499), exp(-6.5185), 0, 0, 0)
initProbs <- initProbs / sum(initProbs)

# Add rows of zeros for individuals enrolled over ten years.
Basis <- rbind( diag(11), matrix(0,5,11))
colnames(Basis) <- c('c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11')

cpts_coefs <- c(-1.242275027,-1.112106935,-0.981938842,-0.851770750,-0.721602657,
		        -0.591434565,-0.461266472,-0.331098380,-0.200930287,-0.070762195, 0.000000000)
deathBias <- Basis %*% cpts_coefs

Basis_amyl <- bs( 50:120, knots=c(55,65,75,90), degree=3, intercept=TRUE)
colnames(Basis_amyl) <- c('c1_amyl','c2_amyl','c3_amyl','c4_amyl','c5_amyl','c6_amyl','c7_amyl','c8_amyl')

cpts_cubic_coefs <- c(-6.580919931,-6.501250885,
 				      -6.427228628,-4.673031945,
			          -3.965086053,-3.111756265,
			          -1.950950409, 3.588549810)
cubicSpline <- Basis_amyl %*% cpts_cubic_coefs


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Collect information about the real data set.
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------


load('synthetic_MCSA_data.rda')
useData <- useData[useData$obstype!=0,]

useData$age <- useData$age + 68
useData$iage <- useData$iage + 68
useData$educ <- useData$educ + 14.03802
useData$ntests <- useData$ntests + 3.153419

ptnum <- unique(useData$ptnum)
N_useData <- length(ptnum)
interObsTime <- NULL
entryAgeDist <- NULL
obsEduc <- NULL
propMale <- 0
propAPOE4 <- 0
propDeaths <- 0
NumObs <- rep(0,N_useData)
for(i in 1:N_useData){   
	subject <- useData[useData$ptnum==ptnum[i],,drop=FALSE]
	
	# The number of observations for each subject.
	NumObs[i] <- nrow(subject)
		
	# The times between observations.
	if(!(2 %in% subject$obstype)){  interObsTime <- c( interObsTime, round( diff(subject$age), 4))  }
	
	# Determine the subject's age at entry.
	entryAgeDist <- c( entryAgeDist, subject$age[1])
	
	# Determine the subject's level of education.
	obsEduc <- c( obsEduc, subject$educ[1])
	
	# Determine whether the subject is male.
	propMale <- propMale + as.integer(subject$male[1]==1)
	
	# Determine whether the subject has an apoe4 allele.
	propAPOE4 <- propAPOE4 + subject$apoe4[1]
	
	# Determine whether the subject's death was observed.
	propDeaths <- propDeaths + as.integer(2 %in% subject$obstype)
}
propMale <- propMale / N_useData
propAPOE4 <- propAPOE4 / N_useData
propDeaths <- propDeaths / N_useData


# Participants with at least one scan.
Scanned_ptnums <- unique(useData$ptnum[ !is.na(useData$lpib) | !is.na(useData$thickness) ])
# Participates with at least one of each scan.
#All_scanned_ptnums <- unique(useData$ptnum[ !is.na(useData$lpib) & !is.na(useData$thickness) ])
propScanned <- length( Scanned_ptnums ) / N_useData
# This is the subset of the study participants who received at least one scan.
Scanned <- useData[ useData$obstype != 2 & useData$ptnum %in% Scanned_ptnums, ]
# Of the scanned subset, this is number of clinical visits.
NumScanned <- nrow(Scanned)

# Construct a table of proportions of observed/not observed scan combinations.
freqTable <- matrix(0,2,2)
colnames(freqTable) <- c('A scan', 'no A scan')
rownames(freqTable) <- c('N scan', 'no N scan')
freqTable[1,1] <- nrow( Scanned[ !is.na(Scanned$lpib) & !is.na(Scanned$thickness), ] ) / NumScanned
freqTable[1,2] <- nrow( Scanned[ is.na(Scanned$lpib) & !is.na(Scanned$thickness), ] ) / NumScanned
freqTable[2,1] <- nrow( Scanned[ !is.na(Scanned$lpib) & is.na(Scanned$thickness), ] ) / NumScanned
freqTable[2,2] <- nrow( Scanned[ is.na(Scanned$lpib) & is.na(Scanned$thickness), ] ) / NumScanned





#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Fill in the states using the transition rate matrix and error matrix estimated on the real data set.
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------



Q <- function(iage,male,educ,apoe4,entryAge){
	
	iyrsEnr = floor( max( iage - entryAge, 0) )
	vals = c(1,iage,male,educ,apoe4)
	
	q1  = exp( vals[3:5] %*% betaMat[1,2:4] + cubicSpline[iage-49] ) # Trans from state 1 to state 2.
	q2  = exp( vals %*% betaMat[2,] )  # Transition from state 1 to state 3.
	q3  = exp( vals %*% betaMat[3,] )  # Transition from state 2 to state 4.
	q4  = exp( vals %*% betaMat[4,] )  # Transition from state 3 to state 4.
	q5  = exp( vals %*% betaMat[5,] )  # Transition from state 3 to state 5.
	q6  = exp( vals %*% betaMat[6,] )  # Transition from state 4 to state 6.
	q7  = exp( vals %*% betaMat[7,] )  # Transition from state 5 to state 6.
	# Note that the time index for the basis spline functions starts at zero.  Hence the +1.
	q8  = exp( vals %*% betaMat[8,] + deathBias[iyrsEnr+1] ) # Trans from state 1 to state 7
	q9  = q8      									  		# Transition from state 2 to state 7.
	q10 = q8   										 		# Transition from state 3 to state 7.
	q11 = q8 										 		# Transition from state 4 to state 7.
	q12 = exp( vals %*% betaMat[12,] ) # Transition from state 5 to state 7.
	q13 = exp( vals %*% betaMat[13,] ) # Transition from state 6 to state 7.

	qmat = matrix(c(  0, q1, q2,  0,  0,  0, q8,
					  0,  0,  0, q3,  0,  0, q9,
					  0,  0,  0, q4, q5,  0,q10,
				  	  0,  0,  0,  0,  0, q6,q11,
				  	  0,  0,  0,  0,  0, q7,q12,
				  	  0,  0,  0,  0,  0,  0,q13,
					  0,  0,  0,  0,  0,  0,  0), ncol=7, byrow=TRUE) 
	diag(qmat) = -rowSums(qmat)
	
	return(qmat)
}


# Set the sample size.  Note that the useData data set has #### subjects.
N <- N_useData
rawData <- NULL
propDeaths_sim <- 0
NumObs_sim <- NULL
for(i in 1:N){
	
	# Sample the gender, as proportional to the useData data set.
	male <- as.integer(runif(1,0,1) < propMale)
	
	# Sample the apoe4, as proportional to the useData data set.
	apoe4 <- as.integer(runif(1,0,1) < propAPOE4)
	
	# Sample for an education level.
	educ <- sample( obsEduc, size=1)
	
	# Sample for an age of enrollment.
	entryAge <- sample( entryAgeDist, size=1)

	# Sample for an initial state, conditional on not demented nor dead by the enrollment age.
	P <- initProbs
	if(entryAge >= 51){
		for(iage in 50:(floor(entryAge)-1)){  P <- P %*% expm( Q(iage,male,educ,apoe4,entryAge) )  }
	}
	P <- P %*% expm( (entryAge - floor(entryAge)) * Q(floor(entryAge),male,educ,apoe4,entryAge) )
	trueState <- sample(1:4, size=1, prob= P[1:4]/sum(P[1:4]) )
	
	# Sample the remainder of the underlying state sequence. 
	age <- entryAge
	time1 <- entryAge
	s <- trueState
	while(s < 7){
		
		# Infinitesimal transition rates.
		qmat <- Q(floor(time1),male,educ,apoe4,entryAge)

		# Possible next states.
		moveToStates <- which(qmat[s,] > 0)

		# Sample the wait times before transition to each of the next possible states.
		waitTimes <- rexp( n=length(moveToStates), rate= qmat[s,moveToStates])

		# If any of the wait times are smaller than dt, then transition to the state with the minimum wait time.
		min_waitTime <- min(waitTimes)
		if(min_waitTime < dt){  s <- moveToStates[ which(waitTimes == min_waitTime) ]  }
		
		time1 <- time1 + dt
		
		age <- c( age, time1)
		trueState <- c( trueState, s)	
	}
	timeOfDeath <- tail(age,1)
	
	# Sample inter-observation times from the useData data set.  Maximum of 12 years in study.
	visitTimes <- NULL
	time2 <- entryAge
	maxObsAge <- entryAge + 12
	prob_remain_in_study <- .9125	 
	randUnif <- 0
	
	while(time2 < min( maxObsAge, timeOfDeath) & randUnif < prob_remain_in_study){
		
		visitTimes <- c( visitTimes, time2)
		time2 <- time2 + sample( interObsTime, size=1)
		randUnif <- runif(1,0,1)	 
	}
	# If death occured before the study ended, then record the time of death.
	if( timeOfDeath < maxObsAge & randUnif < prob_remain_in_study )  visitTimes <- c( visitTimes, timeOfDeath)


	n_i <- length(visitTimes)
	state <- NULL
	for(k in 1:n_i){  state <- c( state, tail( trueState[ age <= visitTimes[k] ], 1))  }
	age <- visitTimes
	  							 
	# Generate MMSE scores using the estimated coeficients from the useData data set.
	mmse <- rep(NA,n_i)
	coef <- c( 1,-0.014554096,-0.089207557, 0.058688151,-0.004627584, 0.057826362)
	intercept <- c( 0.538414931, 0.178435764,-0.672421219,-1.827065032,-2.341332296,-5.208797471) 
	intercept <- intercept -68*coef[2] -14.03802*coef[4] -3.153419*coef[6]
	sigma <- exp(-1.261838927)^.5
	for(k in 1:n_i){  
		ntests <- k
		if(state[k]<7){  
			mu <- c( intercept[state[k]], age[k], male, educ, apoe4, ntests) %*% coef
			mmse[k] <- rnorm( n=1, mean=mu, sd=sigma)  
		}
	}
	
	# Sample whether amyloid and thickness measurements are observed, and if so, then generate the measurements.
	lpib <- rep(NA,n_i)
	thickness <- rep(NA,n_i)
	if( runif(1,0,1) < propScanned ){
		
		for(k in 1:n_i){
		 
			# Sample lpib based on A- states.
			if( state[k] %in% c(1,3,5) ){  
				lpib[k] <- rnorm( n=1, mean=-1.084010116, sd=exp(-2.240498159)^.5)  
			# Sample lpib based on A+ states.
			} else if( state[k] %in% c(2,4,6) ){
				lpib[k] <- rnorm( n=1, mean= 0.110684536, sd=exp(-2.240498159)^.5)  
			}
		 
			# Sample thickness based on N- states.
			if( state[k] %in% c(1,2) ){  
				thickness[k] <- rnorm( n=1, mean=2.655601472, sd=exp(-3.869343130)^.5)  
			# Sample thickness based on N+ states.
			} else if( state[k] %in% c(3,4,5,6) ){
				thickness[k] <- rnorm( n=1, mean=2.498948044, sd=exp(-3.869343130)^.5)  
			}
			
			# Sample whether the lpib and/or thickness measurements are observed.
			obsScan <- sample( 1:4, size=1, prob=freqTable)
			if(obsScan == 1){
				# Do nothing, both scans are observed.
			
			} else if(obsScan == 2){
				thickness[k] <- NA
			
			} else if(obsScan == 3){
				lpib[k] <- NA
			
			} else if(obsScan == 4){
				lpib[k] <- NA
				thickness[k] <- NA
			}	
		}
	} 
	
	
	# Label the true dementia status for each observation.
	DemStatus <- ifelse(state %in% c(1,2,3,4), 1, ifelse(state == 7, NA, 2))

	# Add noise to the dementia status.
	demError <- matrix(c(                1,  exp(-4.712799832),
					   	 exp(-2.163822374),                 1), ncol=2, byrow=TRUE)
	demError <- demError / rowSums(demError)
	
	for(k in 1:n_i){  
		if(!is.na(DemStatus[k])){  DemStatus[k] <- sample( 1:2, size=1, prob=demError[DemStatus[k],])  }
	}	
		
		
	ptnum <- i			
	rawData <- rbind( rawData, data.frame(ptnum,age,male,educ,apoe4,mmse,lpib,thickness,DemStatus,state) )

	if(7 %in% state){  propDeaths_sim <- propDeaths_sim + 1  }
	NumObs_sim <- c( NumObs_sim, n_i)

	#print(i)
}

colnames(rawData) <- c('ptnum','age','male','educ','apoe4','mmse','lpib','thickness','DemStatus','state')
N <- length(unique(rawData$ptnum))
propDeaths_sim <- propDeaths_sim / N



#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Add censored rows.
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------



obstype <- rep(1,nrow(rawData))
obstype[rawData$state == 7] = 2
obstype[!duplicated(rawData$ptnum)] = 3

iage <- floor(rawData$age)


hold <- cbind(rawData,obstype,iage)
hold <- hold[,c('ptnum','age','iage','male','educ','apoe4','mmse','lpib','thickness','DemStatus','obstype')]

tempRow <- rep(0,ncol(hold)+1+11+8)
names(tempRow) <- c('ptnum','age','iage','male','educ','apoe4','mmse','lpib','thickness','DemStatus','obstype',
					'ntests','c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11',
					'c1_amyl','c2_amyl','c3_amyl','c4_amyl','c5_amyl','c6_amyl','c7_amyl','c8_amyl')

num <- 1
demData <- NULL
for(i in unique(rawData$ptnum)){
	
	current <- NULL
	subject <- hold[hold$ptnum==i,,drop=FALSE]
	
    #Determine for each visit how many times a subject has taken the mmse.
   	ntests <- rep(0,nrow(subject))
   	for(k in 1:nrow(subject)){  ntests[k] <- sum( !is.na(subject$mmse[1:k]) )  }
	
	# Years enrolled is an integer which only increases with every full-INTEGER year in the study.
	iyrsEnr <- floor( subject$iage - min(subject$age) )
	iyrsEnr[iyrsEnr<0] <- 0
	
	# Note that the time index for the basis spline functions starts at zero.  Hence the +1.
	subject <- cbind( subject, ntests, Basis[ iyrsEnr+1,, drop=FALSE], Basis_amyl[ subject$iage-49,, drop=FALSE])
	
	#------------------------------------
	censoredAges <- 50:max(subject$age)
	for(t in censoredAges ){

		# If 't' corresponds to an observed age, then the next row will include the observed clinical visit data.
		if(t %in% subject$age){	
			current <- rbind( current, subject[subject$iage==floor(t),]) 
		} else{
		
			# Create a CENSORED row for each subject at each INTEGER year of age.
			tempRow['ptnum'] <- i
			tempRow['age'] <- t
			tempRow['iage'] <- t
			tempRow['obstype'] <- 0
			tempRow['male'] <- subject$male[1]
			tempRow['mmse'] <- NA
			tempRow['educ'] <- subject$educ[1]
			tempRow['apoe4'] <- subject$apoe[1]
			tempRow['lpib'] <- NA
			tempRow['thickness'] <- NA
			tempRow['DemStatus'] <- NA
			
			if(t==50){ tempRow['ntests'] <- 0 } else{ tempRow['ntests'] <- current[nrow(current),'ntests'] }
			
			iyrsEnr <- floor( max( t - min(subject$age), 0) )
		
			tempRow[c('c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11')] <- Basis[ iyrsEnr+1,,drop=FALSE]
			tempRow[c('c1_amyl','c2_amyl','c3_amyl','c4_amyl',
					  'c5_amyl','c6_amyl','c7_amyl','c8_amyl')] <- Basis_amyl[ t-49,, drop=FALSE]
			current <- rbind( current, tempRow)
			
			# If 't' corresponds to an observed INTEGER age, then the subject was observed some time during this age.  According, the next row will include the observed clinical visit data.  Recall that integer age is simply the floor(age).
			if(t %in% subject$iage){ current <- rbind( current, subject[subject$iage==t,]) }
		}

	}
	#------------------------------------
	
	demData <- rbind( demData, current)
	#print(num)
	num <- num+1
}
colnames(demData) <-c('ptnum','age','iage','male','educ','apoe4','mmse','lpib','thickness','DemStatus','obstype',
					   'ntests','c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11',
					   'c1_amyl','c2_amyl','c3_amyl','c4_amyl','c5_amyl','c6_amyl','c7_amyl','c8_amyl')
rownames(demData) <- NULL


meanIage <- round( mean(demData$iage), 0)
# Center using the mean iage for the real data set.
save(meanIage,file=paste0(Dir,'meanIage',args[1],'.rda'))
demData$iage <- demData$iage - 68
demData$age <- demData$age - 68

meanEduc <- mean(demData$educ)
# Center using the mean educ for the real data set.
save(meanEduc,file=paste0(Dir,'meanEduc',args[1],'.rda'))
demData$educ <- demData$educ - 14.03802

meanNtests <- mean(demData[!is.na(demData$mmse),]$ntests)
# Center using the mean ntests for the real data set..
save(meanNtests,file=paste0(Dir,'meanNtests',args[1],'.rda'))
demData$ntests <- demData$ntests - 3.153419 

# Save the sample size.
save(N,file=paste0(Dir,'sampleSize',args[1],'.rda'))

save(demData,file=paste0(Dir,'demData',args[1],'.rda'))



#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Compute the frequencies of observed transitions.
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------



# Transition frequencies for the true data set.
nTrans <- rep(0,13)
for(i in unique(useData$ptnum)){
	
	subject <- useData[useData$ptnum==i,,drop=FALSE]
	
	# Note that DemStatus is NA for observations of death.
	state <- ifelse( !is.na(subject$DemStatus),	
					ifelse( subject$DemStatus==1, 
						   ifelse(subject$lpib < log(1.4 - 1),
						    	  ifelse(subject$thickness > 2.74, 1, 3), 	
								  ifelse(subject$thickness > 2.74, 2, 4) 
						   	   	  ),
						   ifelse(subject$lpib < log(1.4 - 1), 5, 6)
						   ),
					7
					)	
	state[is.na(state)] <- -99

	if(7 %in% state){
		if( max(setdiff(state,7)) == 1){nTrans[8]  = nTrans[8] +1}
		if( max(setdiff(state,7)) == 2){nTrans[9]  = nTrans[9] +1}
		if( max(setdiff(state,7)) == 3){nTrans[10] = nTrans[10] +1}
		if( max(setdiff(state,7)) == 4){nTrans[11] = nTrans[11] +1}
		if( max(setdiff(state,7)) == 5){nTrans[12] = nTrans[12] +1}
		if( max(setdiff(state,7)) == 6){nTrans[13] = nTrans[13] +1}
	}
	
	if(   1 %in% state   &   2 %in% state   ){ nTrans[1] = nTrans[1] +1 } # state 1 to state 2.
	if(   1 %in% state   &   3 %in% state   ){ nTrans[2] = nTrans[2] +1 } # state 1 to state 3.
	if(   2 %in% state   &   4 %in% state   ){ nTrans[3] = nTrans[3] +1 } # state 2 to state 4.
	if(   3 %in% state   &   4 %in% state   ){ nTrans[4] = nTrans[4] +1 } # state 3 to state 4.
	if(   3 %in% state   &   5 %in% state   ){ nTrans[5] = nTrans[5] +1 } # state 3 to state 5.
	if(   4 %in% state   &   6 %in% state   ){ nTrans[6] = nTrans[6] +1 } # state 4 to state 6.
	if(   5 %in% state   &   6 %in% state   ){ nTrans[7] = nTrans[7] +1 } # state 5 to state 6.
}


# Transition frequencies for the simulated data set.
nTrans_sim <- rep(0,13)
obs_demData <- demData[demData$obstype!=0,]
for(i in unique(obs_demData$ptnum)){

	subject <- obs_demData[obs_demData$ptnum==i,,drop=FALSE]
	
	# Note that DemStatus is NA for observations of death.
	state <- ifelse( !is.na(subject$DemStatus),	
					ifelse( subject$DemStatus==1, 
						   ifelse(subject$lpib < log(1.4 - 1),
						    	  ifelse(subject$thickness > 2.74, 1, 3), 	
								  ifelse(subject$thickness > 2.74, 2, 4) 
						   	   	  ),
						   ifelse(subject$lpib < log(1.4 - 1), 5, 6)
						   ),
					7
					)	
	state[is.na(state)] <- -99
	
	if(7 %in% state){
		if( max(setdiff(state,7)) == 1){nTrans_sim[8]  = nTrans_sim[8] +1}
		if( max(setdiff(state,7)) == 2){nTrans_sim[9]  = nTrans_sim[9] +1}
		if( max(setdiff(state,7)) == 3){nTrans_sim[10] = nTrans_sim[10] +1}
		if( max(setdiff(state,7)) == 4){nTrans_sim[11] = nTrans_sim[11] +1}
		if( max(setdiff(state,7)) == 5){nTrans_sim[12] = nTrans_sim[12] +1}
		if( max(setdiff(state,7)) == 6){nTrans_sim[13] = nTrans_sim[13] +1}
	}
	
	if(   1 %in% state   &   2 %in% state   ){ nTrans_sim[1] = nTrans_sim[1] +1 } # state 1 to state 2.
	if(   1 %in% state   &   3 %in% state   ){ nTrans_sim[2] = nTrans_sim[2] +1 } # state 1 to state 3.
	if(   2 %in% state   &   4 %in% state   ){ nTrans_sim[3] = nTrans_sim[3] +1 } # state 2 to state 4.
	if(   3 %in% state   &   4 %in% state   ){ nTrans_sim[4] = nTrans_sim[4] +1 } # state 3 to state 4.
	if(   3 %in% state   &   5 %in% state   ){ nTrans_sim[5] = nTrans_sim[5] +1 } # state 3 to state 5.
	if(   4 %in% state   &   6 %in% state   ){ nTrans_sim[6] = nTrans_sim[6] +1 } # state 4 to state 5.
	if(   5 %in% state   &   6 %in% state   ){ nTrans_sim[7] = nTrans_sim[7] +1 } # state 5 to state 6.
}

obs_demData <- demData[!(demData$obstype %in% c(0,2)),]
Scanned <- obs_demData[ obs_demData$ptnum %in% 
					    unique(obs_demData$ptnum[ !is.na(obs_demData$lpib) | !is.na(obs_demData$thickness) ]),]
# Of the scanned subset, this is number of clinical visits.
NumScanned <- nrow(Scanned)

# Construct a table of proportions of observed/not observed scan combinations.
freqTable_sim <- matrix(0,2,2)
colnames(freqTable_sim) <- c('A scan', 'no A scan')
rownames(freqTable_sim) <- c('N scan', 'no N scan')
freqTable_sim[1,1] <- nrow( Scanned[ !is.na(Scanned$lpib) & !is.na(Scanned$thickness), ] ) / NumScanned
freqTable_sim[1,2] <- nrow( Scanned[ is.na(Scanned$lpib) & !is.na(Scanned$thickness), ] ) / NumScanned
freqTable_sim[2,1] <- nrow( Scanned[ !is.na(Scanned$lpib) & is.na(Scanned$thickness), ] ) / NumScanned
freqTable_sim[2,2] <- nrow( Scanned[ is.na(Scanned$lpib) & is.na(Scanned$thickness), ] ) / NumScanned


cat('useData data set sample size                             = ', N_useData,'\n')
cat('useData data set transition fequencies                   = ', round( nTrans / sum(nTrans), 2),'\n')
cat('useData data set transition counts                       = ', nTrans,'\n')
cat('useData data set proportion of observed deaths           = ', propDeaths,'\n')
cat('useData data set quantiles of number of observations     = ','\n')
print(quantile(NumObs))
cat('useData data set proportion of observed scans = ','\n')
print(freqTable)
cat('\n')
cat('Simulated data set sample size                         = ', N,'\n')
cat('Simulated data set transition fequencies               = ', round( nTrans_sim / sum(nTrans_sim), 2),'\n')
cat('Simulated data set transition counts                   = ', nTrans_sim,'\n')
cat('Simulated data set proportion of observed deaths       = ', propDeaths_sim,'\n')
cat('Simulated data set quantiles of number of observations = ','\n')
print(quantile(NumObs_sim))
cat('Simulated data set proportion of observed scans = ','\n')
print(freqTable_sim)




