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
# https://github.com/emmesgit/RSEQREP/SOFTWARE.xlsx  
# 
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# To cite this software, please reference doi:10.12688/f1000research.13049.1
#
# Program:  init-euclid-dist.r
# Version:  RSEQREP 1.1.0
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Initialize Euclidean distances used for hierachical clustering and 
#				multidimensional scaling
# Input:    analysis/lcpm/<spc>_alltp_lcpm_<norm>_normalized_filtered.tab
# Output:  	analysis/dist/<spc>_alltp_euc_dist_<norm>_normalized_<std>.RData
#############################################################################################################

source('init-analysis.r')

## set directories
in.dir.lcpm  = paste(res.dir,'lcpm',sep='/');
out.dir.dist = paste(res.dir,'dist',sep='/');

## update specimen flag to inclue a variable across all
## types if there is more than 1 specimen type
if (length(spcFlags)>1) {
	spcFlags = c('all',spcFlags)
}

## loop through specimen types
for(s in 1:length(spcFlags)) {
	spcFlag = spcFlags[s];
		
	## loop through normalization types
	for(n  in 1:length(nrmFlags)) {
		nrmFlag=nrmFlags[n];
		
		## import data
		in.file = paste(in.dir.lcpm,'/',spcFlag,'_alltp_lcpm_',nrmFlag,'_normalized_filtered.tab.gz',sep='');
		dta = read.csv(in.file,header=T,stringsAsFactors=F,check.names = FALSE,sep='\t');
		
		## loop through standardization types
		for (d in 1:length(stdFlags)) {
			stdFlag = stdFlags[d];
				
			## align meta data and samples 
			mta.dist = mta[na.omit(match(colnames(dta),mta$samid)),]
			dta.dist = t(dta[,match(mta.dist$samid,colnames(dta))]);
			
			## standardize column variables (mean=0, var=1) if std flag is set
			if(stdFlag=='std') {
				dta.dist=apply(dta.dist,2,scale);
			} 
			
			## calculate eculidean distances;
			dist.res= dist(dta.dist,method='euclidean');
			
			## update rownames to sample ids
			attr(dist.res, "labels") = mta.dist$samid;
			
			## save results as R object
			filename = paste(out.dir.dist,'/',spcFlag,'_alltp_euc_dist_',nrmFlag,'_normalized_',stdFlag,'.RData',sep='');
			save(dist.res,file=filename,compress=T);		
		}
	}
}