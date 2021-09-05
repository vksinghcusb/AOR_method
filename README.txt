##****************************************************************************************####
			 A R-function "modTraitAssociation" 							
##----------------------------------------------------------------------------------------####
Description												
																							 
The R function "modTraitAssociation" quantifies association betweeen co-expression module	 
and clinical trait of interest. The function implement the analysis of rank (AOR) method	
to calulate Gs score and associated p-value for each of the module. The high Gs score 
represent higher expression of module genes in positve class sample whereas low Gs score	 
denote lower expression of module genes in positive class samples. The p-value is obtained 
using random permutaton testing by creating a null distribution of Gs score for concerned	 
module. The function depends on R-package WGCNA. 																					 

Usage:															
 											
modTraitAssociation(datExpr, moduleColors, datTrait, heatmap=FALSE)						 
																				
Arguments:																				 
																							 
datExpr			The expression matrix, sample in rows and genes/probes in columns		 
																							 
moduleColors		A character vector same as length = ncol(ExpData) assigning each		 
			gene/probe to a unique module.											 
																							 
datTraits		A named numeric vector same as length = nrow(Expdata). The values of 	 
			this vector indicate the class to wich the corresponding sample belongs. 
			(Control class samples = 0 and positive class sample =1 ) or a numeric
			vector of length = nrow(ExpData) giving measurements for clinical traits
			of interes for each of the sample.		
																							  
heatmap			If set to TRUE create a heatmap representation of association sores and 
			corresponding p-values obtained by AOR and Eigengene method			
																							 
Value																						 
																							 
A list containg scaled Gs score and corresponsing p-values for each of the module				 
																							 

##############################################################################################
