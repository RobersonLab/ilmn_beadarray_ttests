VERSION="1.0.4"

setClass( "ilmnChipData", representation(exprs="matrix", detect.exprs="matrix", nbeads="matrix", sd.exprs="matrix"))

#########################################################
removeUnexpInTests = function( chipData, colA, colB, pvalCutoff=0.05 )
{
	keepsIndex = as.numeric( sort( unique( c( which(chipData@detect.exprs[,colA] <= pvalCutoff), which(chipData@detect.exprs[,colB] <= pvalCutoff) ) ) ) )
	
	chipData@exprs = chipData@exprs[keepsIndex,]
	chipData@sd.exprs = chipData@sd.exprs[keepsIndex,]
	chipData@detect.exprs = chipData@detect.exprs[keepsIndex,]
	chipData@nbeads = chipData@nbeads[keepsIndex,]
	
	return( chipData )
}

indexesOfUnexpInTests = function( chipData, colA, colB, pvalCutoff=0.05 )
{
	keepsIndex = as.numeric( sort( unique( c( which(chipData@detect.exprs[,colA] <= pvalCutoff), which(chipData@detect.exprs[,colB] <= pvalCutoff) ) ) ) )
	
	return( keepsIndex )
}

keepTheseIndexes = function( chipData, indexesToKeep )
{
	chipData@exprs = chipData@exprs[indexesToKeep,]
	chipData@sd.exprs = chipData@sd.exprs[indexesToKeep,]
	chipData@detect.exprs = chipData@detect.exprs[indexesToKeep,]
	chipData@nbeads = chipData@nbeads[indexesToKeep,]
	
	return( chipData )
}

customQuantileMedianNormalize = function( chipData )
{
	# Pull data
	exp.data = chipData@exprs
	nbeads.data = chipData@nbeads
	sd.data = chipData@sd.exprs
	detect.data = chipData@detect.exprs
	
	# Get dims and build array for resorting
	nrows = dim(exp.data)[1]
	ncols = dim(exp.data)[2]
	
	# set the original order for resorting values later
	orig.order = matrix(nrow=nrows, ncol=ncols, rep(seq(1,nrows),ncols))
	
	# sort the data columns
	for (i in 1:ncols)
	{
		sort.order = order(exp.data[,i])
		
		orig.order[,i] = orig.order[sort.order,i]
		exp.data[,i] = exp.data[sort.order,i]
		nbeads.data[,i] = nbeads.data[sort.order,i]
		sd.data[,i] = sd.data[sort.order,i]
		detect.data[,i] = detect.data[sort.order,i]
	}
	rm(sort.order)
	
	# awful memory copy. easy to implement, not too intensive if data is small. reevaluate if this becomes limiting
	exp.cp = exp.data
	nbeads.cp = nbeads.data
	sd.cp = sd.data
	detect.cp = detect.data
	
	# sort all individual rows
	for (i in 1:nrows)
	{
		sort.order = order(exp.cp[i,])
		
		exp.cp[i,] = exp.cp[i,sort.order]
		nbeads.cp[i,] = nbeads.cp[i,sort.order]
		sd.cp[i,] = sd.cp[i,sort.order]
		detect.cp[i,] = detect.cp[i,sort.order]
	}
	rm(sort.order)
	
	midValue = ceiling( ncols/2 )
	
	q.beads = nbeads.cp[,midValue]
	q.exp = exp.cp[,midValue]
	q.sd = sd.cp[,midValue]
	q.detect = detect.cp[,midValue]
	
	# replace data with normalized data
	for (i in 1:ncols)
	{
		resort.index = orig.order[,i]
		
		chipData@exprs[resort.index,i] = q.exp
		chipData@sd.exprs[resort.index,i] = q.sd
		chipData@nbeads[resort.index,i] = q.beads
		chipData@detect.exprs[resort.index,i] = q.detect
	}
	
	# return
	return(chipData)
}

log2Transform = function( chipData )
{
	chipData@sd.exprs = (1/log(2)) * (chipData@sd.exprs / chipData@exprs)
	chipData@exprs = log2(chipData@exprs)
	 
	return( chipData )
}

shiftData = function (chipData, shiftTo=2.0)
{
	# 2.0 default chosen since log2 is commonly applied.
	for (colIndex in 1:(dim(chipData@exprs)[2]))
	{
		lowest.probe = min(chipData@exprs[,colIndex])
		chipData@exprs[,colIndex] = chipData@exprs[,colIndex] + (shiftTo - lowest.probe)
	}
	
	return( chipData )
}

sortMatrixByRowColNames = function( inputMatrix )
{
	inputMatrix = inputMatrix[,order(colnames(inputMatrix))]
	inputMatrix = inputMatrix[order(rownames(inputMatrix)),]
	return( inputMatrix )
}

getTvaluesFromSummary = function (grp1.val, grp1.variance, grp1.n, grp2.val, grp2.variance, grp2.n)
{
	#Ts = (grp1.val - grp2.val) / sqrt ( (grp1.variance/grp1.beads) + (grp2.variance/grp2.beads) )
	return ( (grp1.val - grp2.val) / sqrt( (grp1.variance/grp1.n) + (grp2.variance/grp2.n) ) )
}

getTvaluesFromLists = function( grp1, grp2 )
{
	return( getTvaluesFromSummary( grp1.val=mean(grp1), grp1.variance=var(grp1), grp1.n=length(grp1), grp2.val=mean(grp2), grp2.variance=var(grp2), grp2.n=length(grp2) ) )
}

unequalVarianceDoF_MoserStevens = function(grp1.n, grp1.variance, grp2.n, grp2.variance)
{	
	# From Moser & Stevens 1992 Homogeneity of variance in the two-sample means test
	# unequal variance t-test degrees of freedom
	# numerator: [(1/n1) + (mu/n2)]^2
	# denominator: [1/(n1^2 * (n1-1))] + [ mu^2 / (n2^2 * (n2-1))]
	# mu: variance2 / variance1
	
	dfs.mu = grp2.variance / grp1.variance
	dfs.numerator = ((1/grp1.n) + (dfs.mu/grp2.n))^2
	dfs.denominator = (1 / ((grp1.n^2) * (grp1.n-1))) + ( (dfs.mu^2) / ((grp2.n^2) * (grp2.n-1)))
	return (dfs.numerator / dfs.denominator)
}

unequalVarianceDoF_WelchSatterthwaite = function(grp1.n, grp1.variance, grp2.n, grp2.variance)
{
	df.numerator = ((grp1.variance/grp1.n)+(grp2.variance/grp2.n))^2
	df.denominator = (((grp1.variance/grp1.n)^2) /(grp1.n-1)) + ( ((grp2.variance/grp2.n)^2)/(grp2.n-1) )
	return( df.numerator / df.denominator )
}

getDoFfromList = function( grp1, grp2, method='moser' )
{
	if (method == 'moser')
	{
		return( unequalVarianceDoF_MoserStevens( grp1.n=length(grp1), grp1.variance=var(grp1), grp2.n=length(grp2), grp2.variance=var(grp2)) )
	} else if (method == 'welch')
	{
		return( unequalVarianceDoF_WelchSatterthwaite( grp1.n=length(grp1), grp1.variance=var(grp1), grp2.n=length(grp2), grp2.variance=var(grp2)) )
	} else
	{
		stop( paste( "Method [", method, "] not recognized", sep="" ) )
	}
}

calculatePvalues = function( Tvalues, DoF )
{
	nom_pvals = pt(Tvalues, DoF, lower.tail=TRUE)
	nom_pvals[which(Tvalues>0)] = pt(Tvalues[which(Tvalues>0)], DoF[which(Tvalues>0)], lower.tail=FALSE)
	nom_pvals = nom_pvals * 2
	nom_pvals[which(nom_pvals > 1.0)] = 1.0
	
	return( nom_pvals )
}

customIlmnTtest = function( raw_object, test, control, removeUnexp=TRUE, log2transform=TRUE, qNormalize=TRUE, qMethod="median", annotate.frame=NULL, showBox=FALSE )
{
	require(beadarray)
	
	# pull obj in and sort
	exp_val = sortMatrixByRowColNames( exprs(raw_object) )
	exp_se = sortMatrixByRowColNames( se.exprs(raw_object) )
	exp_beads = sortMatrixByRowColNames( raw_object@assayData$nObservations )
	exp_pvals = sortMatrixByRowColNames( Detection(raw_object) )
	
	# make data object
	#representation(exprs="matrix", se.exprs="matrix", var.exprs="matrix", detect.exprs="matrix", nbeads="matrix"))
	rawObj = new("ilmnChipData")
	rawObj@exprs = exp_val
	rawObj@sd.exprs = exp_se * sqrt(exp_beads)
	rawObj@detect.exprs = exp_pvals
	rawObj@nbeads = exp_beads
	
	# get the correct columns
	colGrp1 = which(colnames(rawObj@exprs) == test)
	colGrp2 = which(colnames(rawObj@exprs) == control)
	
	# Shift data according to the lowest value to get a minimum of 2.0
	if (min(rawObj@exprs) < 2.0)
	{
		rawObj = shiftData( rawObj, shiftTo=2.0 )
	}
	
	# onlyExpressed function to trim -- if chosen
	if (removeUnexp)
	{
		indexOfExpressedGenes = indexesOfUnexpInTests( rawObj, colGrp1, colGrp2 )
	}
	
	# log2 transform -- if chosen
	if (log2transform)
	{
		rawObj = log2Transform( rawObj )
	}
	
	# normalize
	if (qNormalize)
	{
		if (qMethod=="median")
		{
			norm.obj = customQuantileMedianNormalize( rawObj )
		} else
		{
			stop("qMethod not recognized!")
		}
	} else
	{
		norm.obj = rawObj
	}
	
	# boxplot of expression
	if (showBox)
	{
		split.screen(c(2,1))
		screen(1)
		boxplot(log2(rawObj@exprs), main="Raw Data Distribution")
		
		screen(2)
		if (log2transform)
		{
			boxplot(norm.obj@exprs, main="Normalized Data Distribution")
		} else
		{
			boxplot(log2(norm.obj@exprs), main="Normalized Data Distribution")
		}
	}
	
	# Now actually remove things
	if (removeUnexp)
	{
		norm.obj = keepTheseIndexes( norm.obj, indexOfExpressedGenes )
	}
	
	# get variance
	exp_variance = norm.obj@sd.exprs * norm.obj@sd.exprs
	
	# calculate statistics	
	grp1.val = norm.obj@exprs[,colGrp1]
	grp1.beads = norm.obj@nbeads[,colGrp1]
	grp1.variance = exp_variance[,colGrp1]
	grp1.sd = norm.obj@sd.exprs[,colGrp1]
	
	grp2.val = norm.obj@exprs[,colGrp2]
	grp2.beads = norm.obj@nbeads[,colGrp2]
	grp2.variance = exp_variance[,colGrp2]
	grp2.sd = norm.obj@sd.exprs[,colGrp2]
	
	# calculate T value
	Ts = getTvaluesFromSummary( grp1.val=grp1.val, grp1.variance=grp1.variance, grp1.n=grp1.beads, grp2.val=grp2.val, grp2.variance=grp2.variance, grp2.n=grp2.beads)
	
	# calculate degrees of freedom with unequal variances
	dfs = unequalVarianceDoF_WelchSatterthwaite( grp1.n=grp1.beads, grp1.variance=grp1.variance, grp2.n=grp2.beads, grp2.variance=grp2.variance )
	
	# calculate p-values
	nom_pvals = calculatePvalues( Tvalues=Ts, DoF=dfs )
	
	# correct for multiple tests
	corr_pvals = p.adjust(nom_pvals, method="BH")
	
	# Bonferroni p
	bf_pvals = nom_pvals * length(nom_pvals)
	bf_pvals[which(bf_pvals>1.0)]=1.0
	
	# if log2, transform back
	if (log2transform)
	{
		grp1.val = 2^grp1.val
		grp1.sd = log(2) * grp1.sd * grp1.val
		
		grp2.val = 2^grp2.val		
		grp2.sd = log(2) * grp2.sd * grp2.val
	}
	
	# calculate fold changes
	grp1.se = grp1.sd / sqrt(grp1.beads)
	grp1.ci95 = qt(0.975, df=grp1.beads-1) * grp1.se
	grp2.se = grp2.sd / sqrt(grp2.beads)
	grp2.ci95 = qt(0.975, df=grp2.beads-1) * grp2.se
	
	fc = grp1.val / grp2.val
	fc[which(fc<1)] = -1 / fc[which(fc<1)]
	
	grp1.relError = grp1.ci95 / grp1.val
	grp2.relError = grp2.ci95 / grp2.val
	
	fc.err = fc * sqrt( grp1.relError^2 + grp2.relError^2 )
	
	fc.lower = fc - fc.err
	fc.upper = fc + fc.err
	
	diff_exprs = data.frame(ProbeID=rownames(norm.obj@exprs), test=grp1.val, control=grp2.val, DegreesOfFreedom=dfs, FoldChange=fc, AbsFC=abs(fc), FC_Err=fc.err, FC_Lower=fc.lower, FC_Upper=fc.upper, T=Ts, p.value=nom_pvals, q.value=corr_pvals, Bonferroni=bf_pvals)
	rownames(diff_exprs) = diff_exprs$Probe
	
	if (!is.null(annotate.frame))
	{
		annotate.frame = annotate.frame[,c("ProbeID", "ILMN_GENE")]
		diff_exprs = merge(diff_exprs, annotate.frame, by="ProbeID")
	}
	
	diff_exprs = diff_exprs[order(diff_exprs$q.value, diff_exprs$Bonferroni, diff_exprs$AbsFC),]
	
	if (!is.null(annotate.frame))
	{
		diff_exprs = diff_exprs[, c("ProbeID", "ILMN_GENE", "FoldChange", "q.value", "Bonferroni", "AbsFC", "test", "control", "DegreesOfFreedom", "FC_Err", "FC_Lower", "FC_Upper", "T", "p.value")]
	} else
	{
		diff_exprs = diff_exprs[, c("ProbeID", "FoldChange", "q.value", "Bonferroni", "AbsFC", "test", "control", "DegreesOfFreedom", "FC_Err", "FC_Lower", "FC_Upper", "T", "p.value")]
	}
	
	colnames(diff_exprs)[ which(colnames(diff_exprs) == "test") ] = test
	colnames(diff_exprs)[ which(colnames(diff_exprs) == "control") ] = control
	colnames(diff_exprs)[ which(colnames(diff_exprs) == "FoldChange") ] = paste("FC ", test, "/", control, sep="")
	
	return(diff_exprs)
}
