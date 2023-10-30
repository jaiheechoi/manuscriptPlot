library(SIMICO)
library(ICSKAT)
library(bindata)
library(fastGHQuad)
library(dplyr)
library(tidyr)
library(R.oo)

# special example where one outcome has no effects
# data generation function where only one outcome has zero effects
gen_dat_oe <- function(bhFunInv, obsTimes = 1:4, windowHalf = 0.1, probMiss = .1, n, p, k, tauSq, gMatCausal, xMat, ES,  nKnots) {

  nocol = p + nKnots + 2
  # number of subjects and outcomes
  # true model has nothing
  fixedMat <- matrix(data=0, nrow=n, ncol=k)

  # get genetic vectors
      # calculate the effect size
  B1 <- rep(NA, ncol(gMatCausal))
  for(j in 1:length(B1)){
    MAF <- apply(gMatCausal, 2, function(x) mean(x)/2)

    B1[j] = ES* abs(log10(MAF[j])) * ((-1)^(j + 1))
  }

  B2 <- rep(NA, ncol(gMatCausal))
  for(j in 1:length(B2)){
    MAF <- apply(gMatCausal, 2, function(x) mean(x)/2)

    B2[j] = ES* abs(log10(MAF[j])) * ((-1)^(j + 1))
  }

  B3 <- rep(NA, ncol(gMatCausal))
  for(j in 1:length(B3)){
    MAF <- apply(gMatCausal, 2, function(x) mean(x)/2)

    B3[j] = 0* abs(log10(MAF[j])) * ((-1)^(j + 1))
  }

  geneticVec <- c((gMatCausal %*% B1), (gMatCausal %*% B2), (gMatCausal %*% B3))

  # random unif vector for PIT
  unifVec <- runif(n=n * k)

  # vectorize fixed effects, all first phenotype, then all second phenotype, etc.
  fixedVec <- c(fixedMat)

  # random intercept
  randomInt <- rnorm(n=n, sd = sqrt(tauSq))

  # get full term
  randomIntRep <- rep(randomInt, k)

  # add random effect
  etaVec <- fixedVec + randomIntRep + geneticVec

  # probability integral transform - assumes PH model
  toBeInv <- -log(1 - unifVec) / exp(etaVec)

  # all n*K exact failure times
  exactTimesVec <- bhFunInv(toBeInv)

  # all exact failure times for all K phenotypes
  exactTimesMat <- matrix(data=exactTimesVec, nrow=n, ncol=k, byrow = FALSE)

  # hold left and right intervals data for all K phenotypes
  leftTimesMat <- matrix(data=NA, nrow=n, ncol=k)
  rightTimesMat <- matrix(data=NA, nrow=n, ncol=k)
  obsInd <- matrix(data=NA, nrow=n, ncol=k)
  tposInd <- matrix(data=NA, nrow=n, ncol=k)

  # do visits separately for each phenotype
  nVisits <- length(obsTimes)

  for (pheno_it in 1:k) {
            # 1 - probMiss is the chance of making it to the visit
  nVisits <- length(obsTimes)
  madeVisit <- matrix(data = stats::rbinom(n=n*nVisits, size=1, prob=(1 - probMiss)), nrow=n, ncol=nVisits)

  # make sure there is at least one visit for each subject
  nMadeVisits <- apply(madeVisit, 1, sum)
  zeroVisits <- which(nMadeVisits == 0)

  while (length(zeroVisits) > 0) {
    madeVisit[zeroVisits, ] <- matrix(data = stats::rbinom(n=length(zeroVisits) * nVisits, size=1,
                prob=(1 - probMiss)), nrow=length(zeroVisits), ncol=nVisits)

    nMadeVisits <- apply(madeVisit, 1, sum)
    zeroVisits <- which(nMadeVisits == 0)
  }

    # actual visit time is uniformly distributed around the intended obsTime, windowHalf on each side
    visitTime <- sweep(matrix(data = stats::runif(n=n*nVisits, min=-windowHalf, max=windowHalf), nrow=n, ncol=nVisits), MARGIN=2, STATS=obsTimes, FUN="+")

    # get all visits for each subject
    allVisits <- madeVisit * visitTime

    # make the interval for each subject
    allInts <- t(mapply(FUN=ICSKAT::createInt, obsTimes = data.frame(t(allVisits)), eventTime=exactTimesMat[, pheno_it]))

    leftTimesMat[, pheno_it] <- allInts[, 1]
    rightTimesMat[, pheno_it] <- allInts[, 2]

    # event time indicators
    obsInd[, pheno_it] <- ifelse(rightTimesMat[, pheno_it] == Inf, 0, 1)
    tposInd[, pheno_it] <- ifelse(leftTimesMat[, pheno_it] == 0, 0, 1)
  }

 leftArray <- array(data=NA, dim=c(n, p + nKnots + 2, k))

  rightArray <- array(data=NA, dim=c(n, p + nKnots + 2, k))


  for (pheno_it in 1:k) {

    tempDmats <- ICSKAT::make_IC_dmat(xMat = xMat, lt = leftTimesMat[, pheno_it],

                                      rt = rightTimesMat[, pheno_it], obs_ind = obsInd[, pheno_it],

                                      tpos_ind = tposInd[, pheno_it], nKnots=nKnots)
    leftArray[, , pheno_it] <- tempDmats$left_dmat
    rightArray[, , pheno_it] <- tempDmats$right_dmat

  }

  #make placeholder to change the data later
  # *** this is how the functions take the inputs
  dataph <- matrix(NA, nrow= n, ncol = nocol)
  vecN <- rep(NA, n)
  dmatph <- list(right_dmat = dataph, left_dmat = dataph)
  #xmatph <- list(dmats = dmatph, lt = vecN, rt = vecN)
  xmatph <- list(dmats = dmatph)
  allph <- matrix(NA, nrow = n, ncol = k)
  threedmat <- list()
  for(i in 1:k){
    threedmat[[i]] <- xmatph
  }
  samp <- list(xDats = threedmat, ts_all = allph, ob_all = allph)

  for(pheno in 1:k){
    samp$xDats[[pheno]]$dmats$right_dmat <- rightArray[,,pheno]
    samp$xDats[[pheno]]$dmats$left_dmat <- leftArray[,,pheno]
    #samp$xDats[[pheno]]$lt <- leftTimesMat[,pheno]
    #samp$xDats[[pheno]]$rt <- rightTimesMat[,pheno]
  }

  samp$ob_all <- obsInd
  samp$ts_all <- tposInd


  # return
  return(list(exactTimesMat = exactTimesMat, leftTimesMat = leftTimesMat,
              rightTimesMat = rightTimesMat, obsInd = obsInd, tposInd = tposInd, fullDat = samp))

}


# function to attain power simulation
# power function
# inputs are:
# - runs: number of iterations
# - n: total sample size
# - q: total number of genetic variants in SNP set
# - rho: correlation of genetic variant
# - ncausal: number of causal genetic variants
# - k: total number of outcomes
# - effectSizes: k-length vector of genetic effect
# - numberKnots: number of internal knots for spline estimation
testPower <- function(runs, n, q, rho, ncausal, k, effectsizes, numberKnots){

multP <- rep(NA, runs) # multiple outcome (Q_ind)
multSig <- rep(NA, runs)
multQ <- rep(NA, runs) 
multburdQ <- rep(NA,runs) # multiple burden test (Q_cor)
multburdP <- rep(NA,runs)
multburdSig <- rep(NA, runs) 
omniP <- rep(NA,runs) # multiple omnibus test (Q_omni)
omniSig <- rep(NA, runs)
singP <- rep(NA, runs) # single ICSKAT bonferroni
singSig <- rep(NA, runs)
burdenP <- rep(NA, runs) # single ICBurden bonferroni
burdenSig <- rep(NA, runs)
singP_ACAT <- rep(NA, runs) # single ICSKAT ACAT
singACAT_Sig <- rep(NA, runs)
burdenP_ACAT <- rep(NA, runs) # single ICBurden ACAT
burdenACAT_Sig <- rep(NA, runs)
indivICSKAT <- matrix(NA, nrow = runs, ncol = k) # matrix of single ICSKAT test results
indivICBurd <- matrix(NA, nrow = runs, ncol = k) # matrix of single IC burdent test results
lcrAll <- matrix(NA, nrow = runs, ncol = k) # left censoring rate for each outcome
rcrAll <- matrix(NA, nrow = runs, ncol = k) # right censoring rate for each outcome
icrAll <- matrix(NA, nrow = runs, ncol = k) # interval censoring rate for each outcome
didConverge <- rep(1, runs) # indicator for whether newton-raphson converged


for(i in 1:runs){
	# Fixed effects
	xMat <- cbind(rnorm(n, mean = 0, sd = 2), rbinom(n=n, size=1, prob=0.5))
	p = ncol(xMat)

	# Genetic effects
	gMat <- sim_gMat(n, q, rho, maxMP = 0.05)

	# Get indices to specific select causal variants
	idx <- Get_CausalSNPs_bynum(gMat, ncausal, 0.05)

	# Subset the gMat
	gMatCausal <- gMat[,idx]

	# True model has nothing
	fixedMat <- matrix(data=0, nrow=n, ncol=k)
        
	# for real times
	bhFunInv <- function(x) {x}
	
	# Generate the multiple outcomes
	exampleDat <- gen_dat_oe(bhFunInv = bhFunInv, obsTimes = 1:4,
			windowHalf = 0.1, probMiss = 0.1, n = n, p = p, k = k, 
	                tauSq = 1, gMatCausal = gMatCausal,
			xMat = xMat, effectSizes = ES, oppSign = TRUE, nKnots = numberKnots)

	# get left and right censoring rates
        lcrAll[i,] <- 1 - colMeans(exampleDat$tposInd)
        rcrAll[i,] <- 1 - colMeans(exampleDat$obsInd)

        for(outs in 1:k){
                leftCtemp <- exampleDat$tposInd[,outs]
                rightCtemp <- exampleDat$obsInd[,outs]
                intCtemp <- rep(0, length(leftCtemp))
                intCtemp[which(leftCtemp == 1 & rightCtemp == 1)] <- 1
                icrAll[i, outs] <- mean(intCtemp)
       }


	# Set the initial estimate values
        nocol = numberKnots + p + 2
	temp_init <- rep(0, nocol)
        temp_init[length(temp_init) - numberKnots] <- 1
        init_beta <- c(rep(temp_init, k),1)

	# Run the newton-raphson
	skip_to_next <- FALSE
	nullFit <- NA
	tryCatch(nullFit <- simico_fit_null(init_beta = init_beta, epsilon = 10^-5, 
                           xDats = exampleDat$fullDat$xDats, 
                           lt_all = exampleDat$leftTimesMat, rt_all = exampleDat$rightTimesMat,
                           k = k, d = 100), error = function(e){skip_to_next <<- TRUE})

        if(is.na(nullFit[1])) {didConverge[i] <- 0}
        if(skip_to_next){ next}

	# Get the test statistics p-values
	out <- simico_out(nullFit = nullFit$beta_fit, xDats = exampleDat$fullDat$xDats, 
                  lt_all = exampleDat$leftTimesMat, rt_all = exampleDat$rightTimesMat, 
                  Itt = nullFit$jmat, a1 = 1, a2 = 25, 
                  G = gMat, k  = k, d = 100)	
	
	#now test individual outcomes
	single_pvals <- rep(NA, k)
	burden_pvals <- rep(NA, k)
	adjusted_pvals <- rep(NA,k)
	wu_pvals <- rep(NA,k)
	wu_pvals_trunc <- rep(NA, k)
	# loop through all the outcomes
    	for(m in 1:k) 
	{
		dmats <- exampleDat$fullDat$xDats[[m]]$dmats
		lt <- exampleDat$leftTimesMat[,m]
      		rt <- exampleDat$rightTimesMat[,m]
      		obs_ind <- exampleDat$fullDat$ob_all[,m]
      		tpos_ind <- exampleDat$fullDat$ts_all[,m]
    		nocol <- ncol(dmats$left_dmat)
  
		# fit null model
      		nullFit <- ICSKAT_fit_null(init_beta = rep(0, nocol), left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat, obs_ind = obs_ind, tpos_ind = tpos_ind, lt = lt, rt = rt)

      		# perform the ICSKAT and Burden tests
     		if(is.na(nullFit$beta_fit)){
     			icskatOut <- list(p_SKAT = NA, p_burden = NA)} else{
      			icskatOut <- ICskat(left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat, lt = lt, rt = rt, obs_ind = obs_ind, tpos_ind = tpos_ind, gMat = gMat, null_beta = as.numeric(nullFit$beta_fit), Itt = nullFit$Itt) }

			#save the pvals
      			single_pvals[m] <- icskatOut$p_SKAT
                	burden_pvals[m] <- icskatOut$p_burden


	} #end of looping though k
	alpha = .05
	
	# bonferroni
	singP[i] <- min(single_pvals)
	singSig[i] <- as.numeric(min(single_pvals) < alpha/k)

	# multiple outcome test
    	multQ[i] <- out$multQ
    	multP[i] <- out$multP
        multSig[i] <- as.numeric(out$multP < alpha)
	# multiple burden test
	multburdQ[i] <- out$burdQ
    	multburdP[i] <- out$burdP
	multburdSig[i] <- as.numeric(out$burdP < alpha)
    
    	# multiple omnibus
	omniP[i] <- out$omniP
	omniSig[i] <- as.numeric(out$omniP < alpha)

	# individual burden bonferroni
	burdenP[i] <- min(burden_pvals)
    	burdenSig[i] <- as.numeric(min(burden_pvals) < alpha/k)

	if(is.na(single_pvals[1])){
		singP_ACAT[i] <- NA
    		singACAT_Sig[i] <- NA
    		burdenP_ACAT[i] <- NA
   		burdenACAT_Sig[i] <- NA 
	} else {
    		singP_ACAT[i] <- ACAT(single_pvals)
		singACAT_Sig[i] <- as.numeric(singP_ACAT[i] < alpha)
    		burdenP_ACAT[i] <- ACAT(burden_pvals)
		burdenACAT_Sig[i] <- as.numeric(burdenP_ACAT[i] < alpha)
		indivICSKAT[i,] <- single_pvals
		indivICBurd[i,] <- burden_pvals
    		}

  } # end of iterations



 sim_results <- data.frame(cbind(multQ, multP, multSig, multburdQ, multburdP, multburdSig, omniP, omniSig, singP, singSig, burdenP, burdenSig, singP_ACAT, singACAT_Sig, burdenP_ACAT, burdenACAT_Sig, indivICSKAT, indivICBurd, lcrAll, rcrAll, icrAll, didConverge))

icrIdx2 <- ncol(sim_results) - 1
icrIdx1 <- icrIdx2 - k + 1
rcrIdx2 <- icrIdx1 - 1
rcrIdx1 <- rcrIdx2 - k + 1
lcrIdx2 <- rcrIdx1 - 1
lcrIdx1 <- lcrIdx2 - k + 1
burIdx2 <- lcrIdx1 - 1
burIdx1 <- burIdx2 - k + 1
indIdx2 <- burIdx1 - 1
indIdx1 <- indIdx2 - k + 1

colnames(sim_results)[indIdx1:indIdx2] <- paste("IndividICSKAT", 1:k, sep = "")
colnames(sim_results)[burIdx1:burIdx2] <- paste("IndivBurden", 1:k, sep = "")
colnames(sim_results)[lcrIdx1:lcrIdx2] <- paste("lcr", 1:k, sep = "")
colnames(sim_results)[rcrIdx1:rcrIdx2] <- paste("rcr", 1:k, sep = "")
colnames(sim_results)[icrIdx1:icrIdx2] <- paste("icr", 1:k, sep = "")


#colnames(sim_results)[15:11] <- paste("IndividICSKAT", 1:3, sep = "")
#colnames(sim_results)[12:14] <- paste("IndivBurden", 1:3, sep = "")
  return(list(sim_results =sim_results, outcomes = k, numCV = ncausal, ES = effectsizes, nKnots = numberKnots))
} # end of function


############################################################################
############################################################################
############################################################################
runs = 200
n = 2400
q = 50 
rho = 0.1
causalVariants = 4
k = 3
numberKnots = c(1)
esVec <- c(0.05, .1, .15, .2, .25, .3, .35, .4)

mat <- NA

# run the test for 200 iterations
for(r in 1:length(esVec)){
  ES <- rep(esVec[r], k)
  temp <- testPower(runs, n, q, rho, MAFcut, causalVariants, k, ES, numberKnots)
  temp_withKnotNum <- cbind(temp$sim_results, "num" = causalVariants, "ES" = ES[1], n)
  mat <- rbind(mat, temp_withKnotNum)
}

mat <- mat[-1,]
