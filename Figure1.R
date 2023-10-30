library(SIMICO)
library(ICSKAT)
library(bindata)
library(fastGHQuad)
library(dplyr)
library(tidyr)
library(R.oo)

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
	exampleDat <- simico_gen_dat(bhFunInv = bhFunInv, obsTimes = 1:4,
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
n = 24000
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
