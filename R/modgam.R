modgam <-
function (rdata, rgrid, permute=0, m="adjusted", sp=NULL, keep=F, verbose=T, ...) {
	x=length(rdata)
	if (x < 3) stop("rdata must include at least 3 columns")
	if (tolower(m)=='adjusted') {
		m = "Adjusted"
		if (x==3) stop("cannot run adjusted model without additional covariates in rdata")
		# fill in median values (or first mode for factors) for any missing covariates in rgrid 
		addvars = names(rdata)[!names(rdata) %in% names(rgrid)][-1]
		if (length(addvars)>0) {
			for (i in addvars) {
				if (is.factor(rdata[,i]))  rgrid[i]=factor(names(which.max(table(rdata[,i]))),levels=levels(rdata[,i])) else
					rgrid[i]=quantile(rdata[,i],0.5)		
				cat(paste("GAM predictions will use ",i," = ",rgrid[1,i]," at all grid points.",sep=""),fill=T) 
			}
		}
	} else
	if (tolower(m)=='unadjusted' | tolower(m)=='crude') {
		x = 3
		m = "Unadjusted"
	} else
		stop(paste("model type",m,"not recognized"))
	if (is.null(sp)) sp=optspan(rdata,m)
	# Creates the crude or adjusted model using the names of variables in the data set, rgrid. Also creates the reduced model
	#	omitting the smooth term
	if (x == 3) {fmla = as.formula(paste(names(rdata)[1], paste("lo(",paste(names(rdata)[2:3], collapse=","),",span=",sp,")"), sep="~"));
	fmla.0 = as.formula(paste(names(rdata)[1], paste("1"), sep="~"))} else 
	{fmla = as.formula(paste(names(rdata)[1], paste("lo(",paste(names(rdata)[2:3], collapse=","),",span=",sp,")+",paste(names(rdata)[-(1:3)], collapse= "+")), sep="~"));
	fmla.0 = as.formula(paste(names(rdata)[1], paste(names(rdata)[-(1:3)], collapse= "+"), sep="~"))}
	if (verbose) {
		cat(paste("The ",tolower(m)," model is: ",sep=""),fill=T) 
		print(fmla,showEnv=F)
	} 
	model=gam(fmla, family = binomial(logit), data = rdata, ...) 
	model.0=gam(fmla.0, family = binomial(logit), data = rdata, ...)
	nullmod=ifelse (x == 3, mean(predict.gam(model.0)), mean(predict.gam(model.0,rgrid))) # aspatial model for calculating ORs

	# Calculates the predicted log odds
	origresults = predict.gam(model,rgrid)					# needs rgrid to include covariates 
	# Converts from log odds to odds ratios (ORs) using the whole study area as the reference, 
	# dividing the odds at each grid point by the odds calculated by the reduced model omitting the smoothing term.
	# store grid, ORs, and call arguments in results list  
	results=list(grid=rgrid[,1:2],OR=as.vector(exp(origresults-nullmod)),m=m,span=sp,gamobj=model)

	# Runs permutation tests
	if (permute>0) { 	##PERMUTATION TEST##
		n=length(rgrid[,1])
		nobs = length(rdata[,1])
		permresults=NULL 						# for permuted OR estimates
		ptranks=rep(1, n)						# for pointwise ranks
		devstat=rep(NA, permute)						# for global test
		devstat[1]=anova(model.0, model)$Deviance[2]	# deviance statistic for original data
		coords=rdata[,2:3]  			
		m.data=rdata
		for (i in 2:permute) {
			index=sample(1:nobs, replace=F)
			m.data[,2:3]=coords[index,]				# randomly reassign individuals to the eligible residences
			# For each permutation, we run the GAM (using a fixed span size)
			m.gam=gam(fmla, family = binomial, data = m.data, ...)
			devstat[i]=anova(model.0, m.gam)$Deviance[2]	# deviance statistic for the permuted data
			tempresults=predict.gam(m.gam, rgrid)		# estimate log odds at each point from the permuted data
			ptranks = ptranks + (origresults > tempresults) 
			if (keep) permresults=cbind(permresults,tempresults)	# add log odds estimate vector to stored results		 
			if (verbose && i%%10==0) cat(paste("Permutation",i,"of",permute),fill=TRUE)
			}
		devglobp=(permute-rank(devstat)[1])/permute			# global p-value based on permuted deviance statistic
        globprint = if (devglobp == 0) paste("<",round(1/permute,3),sep="") else devglobp
		cat(paste("The global statistic for the ",tolower(m)," model is ",globprint,sep=""), fill=TRUE)
		results$global = devglobp								# save global percentile rank
		results$pointwise = as.vector(ptranks/permute)			# save pointwise percentile ranks
		if (keep) {
			results$permutations=exp(permresults-nullmod)		# save permuted odds ratios
			results$globaldevs=devstat							# save permuted global deviance statistics
			}
		}
	return(results)
	}