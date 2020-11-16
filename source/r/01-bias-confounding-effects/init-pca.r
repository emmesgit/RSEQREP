#############################################################################################################
# RSEQREP: RNA-Seq Reports, an open-source cloud-enabled framework for reproducible
# RNA-Seq data processing, analysis, and result reporting
# 
# https://github.com/emmesgit/RSEQREP
#
# Copyright (C) 2017 The Emmes Corporation 
# 
# This program is free software that contains third party software subject to various licenses, 
# namely, the GNU General Public License version 3 (or later), the GNU Affero General Public License 
# version 3 (or later), and the LaTeX Project Public License v.1.3(c). A list of the software contained 
# in this program, including the applicable licenses, can be accessed here: 
# 
# https://github.com/emmesgit/RSEQREP/blob/master/SOFTWARE.xlsx  
# 
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# To cite this software, please reference doi:10.12688/f1000research.13049.1
#
# Program:  init-pca.r
# Version:  RSEQREP 2.0.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Run PCA (gene variables) for the first three components; 
#				run results for normalized and unnormalized, 
#				and standardized and unstandardized data.
# Input:    analysis/lcpm/<spc>_alltp_lcpm_<norm>_normalized_filtered.tab
# Output:  	analysis/pca/<spc>_<time>_pca_<norm>_normalized_<std>.RData
#############################################################################################################

source('init-analysis.r')

## set directories
in.dir.lcpm  = paste(res.dir,'lcpm',sep='/');
out.dir.pca = paste(res.dir,'pca',sep='/');

## update specimen flag to inclue a variable across all
## types if there is more than 1 specimen type
if (length(spcFlags)>1) {
	spcFlags = c('all',spcFlags)
}

## loop through specimen types
for(s in 1:length(spcFlags)) {
	spcFlag = spcFlags[s];
	
	## loop through normalization types
	for(n in 1:length(nrmFlags)) {
		nrmFlag = nrmFlags[n];
		
		## read data
		in.file = paste(in.dir.lcpm,'/',spcFlag,'_alltp_lcpm_',nrmFlag,'_normalized_filtered.tab.gz',sep='');
		dta = read.csv(in.file,header=T,stringsAsFactors=F,check.names = FALSE,sep='\t');
		
		## loop through standardization types
		for(d in 1:length(stdFlags)) {
			stdFlag = stdFlags[d];
		
			## align meta data and samples 
			mta.pca=mta[na.omit(match(colnames(dta),mta$samid)),]
		
			dta.pca = t(dta[,match(mta.pca$samid,colnames(dta))]);
			
					
			## standardize column variables (mean=0, var=1) if std flag is set
			if(stdFlag=='std') {
				dta.pca=apply(dta.pca,2,scale);
				
				## remove rows with NA throught
				not.na.idx = which(apply(dta.pca,2,function(x){!any(is.na(x))}))
				dta.pca = dta.pca[,not.na.idx]
				dta = dta[not.na.idx,]
			} 
			
			## update rownames to sample ids
			## and colnames to gene ids
			rownames(dta.pca)=mta.pca$samid;
			colnames(dta.pca)=rownames(dta);
			
			## execute pca on gene variables
			pca.res = prcomp(dta.pca,scale=F);
			
			## store only the first three components (coefficents and ratated values)
			pca.res$rotation = pca.res$rotation[,1:3];
			pca.res$x = pca.res$x[,1:3];
					
			## calculate percent explained variance
			pca.res$pvar=round(pca.res$sdev^2/sum(pca.res$sdev^2)*100,1);
				
			## save results as R object
			filename = paste(out.dir.pca,'/',spcFlag,'_alltp_pca_',nrmFlag,'_normalized_',stdFlag,'.RData',sep='');
			save(pca.res,file=filename,compress=T);
			
		}
	}
}