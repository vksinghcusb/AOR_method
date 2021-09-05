##****************************************************************************************####
## A R-function "modTraitAssociation" 									              	  ####
##----------------------------------------------------------------------------------------####
# Description																				 #
#																							 #
# The R function "modTraitAssociation" quantifies association betweeen co-expression module	 #
# and clinical trait of interest. The function implement the analysis of rank (AOR) method	 #
# to calulate Gs score and associated p-value for each of the module. The high Gs score 	 #
# represent higher expression of module genes in positve class sample whereas low Gs score	 #
# denote lower expression of module genes in positive class samples. The p-value is obtained #
# using random permutaton testing by creating a null distribution of Gs score for concerned	 #
# module. The function depends on R-package WGCNA. 																					 #
#																							 #
# Usage:																					 #
#																							 # 											
# modTraitAssociation(datExpr, moduleColors, datTrait, heatmap=FALSE)						 #
#																							 #
# Arguments:																				 #
#																							 #
# datExpr			The expression matrix, sample in rows and genes/probes in columns		 #
#																							 #
# moduleColors		A character vector same as length = ncol(ExpData) assigning each		 #
#					gene/probe to a unique module.											 #
#																							 #
# datTraits			A named numeric vector same as length = nrow(Expdata). The values of 	 #
#					this vector indicate the class to wich the corresponding sample belongs. #
#					(Control class samples = 0 and positive class sample =1 )				 #
#					or a numeric vector of length = nrow(ExpData) giving measurements for 	 #
#					clinical traits of interes for each of the sample.		 				 #
#																							 # 
# heatmap			If set to TRUE create a heatmap representation of association sores and  #
#					corresponding p-values obtained by AOR and Eigengene method				 #
#																							 #
# Value																						 #
#																							 #
# A list containg scaled Gs score and corresponsing p-values for each of the module				 #
#																							 #
##############################################################################################

#***************************************************************#
# An R function to reorder the indices of the position vector	#	 
# This function is used by main function "modTraitAssociation"	#
#																#
orderOFindices <- function(x){									#
	y <- sort(unique(x),decreasing=T);							#
																#			
	pVec <- numeric();											#
																#
	for(i in 1:length(y)){										#
		z <- which(x==y[i]);										#
		L <- length(z);											#
		if(L < 3){												#
			pVec <- c(pVec,z);									#
		} else {												#
			iVec <- numeric(length=L);							#
			iVec[seq(1,L,2)] <- 1:length(seq(1,L,2));			#
			iVec[seq(2,L,2)] <- L:(length(seq(1,L,2))+1);		#
			pVec <- c(pVec,z[iVec]);							#
		}														#
																#		
	}															#
return(pVec);													#
}																#
#################################################################






#***************************************************************#
# An R function to arrange sample according to expression value #
# of a given gene												#	 
# This function is used by main function "modTraitAssociation"	#

SamplOrder <- function(vec,datTraits){							#
	return(datTraits[names(sort(vec))])							#
}																#
#################################################################





#***************************************************************************************#
# An R function to calculate Gs score and associated p-values							#
#																						#
modTraitAssociation <- function(datExpr, moduleColors, datTraits, heatmap=FALSE){
library(WGCNA)

	if(!class(datTraits)%in%c("numeric","integer")){
		stop(paste('clinical trait is not a numeric vector'))
	}
	if(is.null(names(datTraits))){
		stop(paste('clinical trait vector is not a named vector'))
	}
	t_1 <- sum(sort(unique(datTraits))==c(0,1))
		
	uniqModules <- unique(moduleColors)

	if (t_1 != 2){
		datTrai0 <- as.numeric(datTraits >= median(datTraits))
	} else {
		datTrait0 <- datTraits
	}			 

		numPositSample <- sum(datTraits0)
		sampleSize <- length(datTraits0)
		

		##### Create NUl Gs score distibution for each of the module
		nulGs <- list()	
		for(module in uniqModules){	
			inx <- which(moduleColors!=module)	
			Mlen <- sum(moduleColors==module)
			tempGs <- numeric()	
			for(i in 1:1000){
 				matx <- t(datExpr[,sample(x=inx,size=Mlen,replace=FALSE)])
		 		iMatx <- t(apply(X=matx,MARGIN=1,FUN=SamplOrder,datTraits=datTraits0));
				pVec1 <- orderOFindices(colSums(iMatx));
				tempGs <- c(tempGs,sum(pVec1[1:numPositSample]))
			}
		nulGs[[module]] <- tempGs;		# the list nulGs contains the permutation testing Gs score for each module
		}

	
		
	## Calculate Gs score for modules
	Gs <- numeric()	
	for(module in uniqModules){
		inx <- which(moduleColors==module);
		matx <- t(datExpr[,inx]);
		iMatx <- t(apply(X=matx,MARGIN=1,FUN=SamplOrder,datTraits=datTraits0)); 
		pVec1 <- orderOFindices(colSums(iMatx));
		Gs <- c(Gs,sum(pVec1[1:numPositSample]))
	}

	## Calculate P-values associated with module Gs score
	nulGsMean <- unlist(lapply(X=nulGs,FUN=mean))
	nulGsSD <- unlist(lapply(X=nulGs,FUN=sd))
	pVal <- numeric()
	for(i in 1:length(uniqModules)){
		if(Gs[i]<=nulGsMean[i]){
			pVal <- c(pVal,pnorm(q=Gs[i],mean=nulGsMean[i],sd=nulGsSD[i],lower.tail=TRUE)*2)
		} else {
			pVal <- c(pVal,pnorm(q=Gs[i],mean=nulGsMean[i],sd=nulGsSD[i],lower.tail=FALSE)*2)
		}
	}


## Quantifying module{trait associations using Eigengene method
	# Define numbers of genes and samples
	nGenes = ncol(datExpr);
	nSamples = nrow(datExpr);
	# Recalculate MEs with color labels
	MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
	MEs = orderMEs(MEs0)
	names(MEs) <- sub(pattern="..","",names(MEs))
	moduleTraitCor = cor(MEs, datTraits0, use = "p");
	moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
	
	Gs_max <- sum(sampleSize:(sampleSize-numPositSample+1))
	Gs_min <- sum(1:numPositSample)
	Gs_scaled <- 2*((Gs-Gs_min)/(Gs_max-Gs_min)) -1

	moduleTraitCor <- as.matrix(cbind(moduleTraitCor[uniqModules,1],(Gs-nulGsMean)/nulGsSD,Gs_scaled))
	moduleTraitPvalue <- as.matrix(cbind(moduleTraitPvalue[uniqModules,1],pVal))
	
	## graphical representation of module trait association
	if(heatmap == TRUE){

# jpeg("ModuleTraiCoorelation.jpeg", width=6, height=8, unit="in", res=600)
		sizeGrWindow(5.5,8)
		# Will display correlations and their p-values
		textMatrix = paste(signif(moduleTraitCor[,1], 2), "\n(",
		signif(moduleTraitPvalue[,1], 1), ")", sep = "");
#		dim(textMatrix) = dim(moduleTraitCor[,1])
		dim(textMatrix) = c(length(uniqModules),1)
		par(mar = c(4, 6.5, 1, 2));
		par(mfcol=c(1,3))
		# Display the correlation values within a heatmap plot
		labeledHeatmap(Matrix = as.matrix(moduleTraitCor[,1]),
					xLabels = c("\n\n\n\nEigengene method"),
					xLabelsPosition = "bottom",
					xLabelsAngle = 0, xLabelsAdj = 0.5,
					yLabels = uniqModules,
					ySymbols = uniqModules,
					colorLabels = FALSE,
					colors = greenWhiteRed(50),
					textMatrix = textMatrix,
					setStdMargins = FALSE,
					cex.text = 0.4,
					cex.lab = 0.7,
					zlim = c(-1,1)
					)
			

		textMatrix = paste(signif(moduleTraitCor[,2], 2), "\n(",
		signif(moduleTraitPvalue[,2], 1), ")", sep = "");
#		dim(textMatrix) = dim(moduleTraitCor[,1])
		dim(textMatrix) = c(length(uniqModules),1)
		# Display the correlation values within a heatmap plot
		labeledHeatmap(Matrix = as.matrix(moduleTraitCor[,2]),
					xLabels = c("\n\n\n\nAOR method\nGs_standersised"),
					xLabelsPosition = "bottom",
					xLabelsAngle = 0, xLabelsAdj = 0.5,					
					yLabels = uniqModules,
					ySymbols = uniqModules,
					colorLabels = FALSE,
					colors = greenWhiteRed(50),
					textMatrix = textMatrix,
					setStdMargins = FALSE,
					cex.text = 0.4,
					cex.lab = 0.7,
					zlim = c(-3.2,3.2)
					)
					
		textMatrix = paste(signif(moduleTraitCor[,3], 2));   #, "\n(",
#		signif(moduleTraitPvalue[,2], 1), ")", sep = "");
#		dim(textMatrix) = dim(moduleTraitCor[,1])
		dim(textMatrix) = c(length(uniqModules),1)
		# Display the correlation values within a heatmap plot
		labeledHeatmap(Matrix = as.matrix(moduleTraitCor[,3]),
					xLabels = c("\n\n\n\nAOR method\nGs_scaled"),
					xLabelsPosition = "bottom",
					xLabelsAngle = 0, xLabelsAdj = 0.5,
					yLabels = uniqModules,
					ySymbols = uniqModules,
					colorLabels = FALSE,
					colors = greenWhiteRed(50),
					textMatrix = textMatrix,
					setStdMargins = FALSE,
					cex.text = 0.4,
					cex.lab = 0.7,
					zlim = c(-1,1)
					)

#	 dev.off()
	} else {
		return(cbind(uniqModules,Gs_scaled,pVal));
	}

}



